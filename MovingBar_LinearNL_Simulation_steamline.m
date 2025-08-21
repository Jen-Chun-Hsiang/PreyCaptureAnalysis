%% Modified LN Model Implementation for Memory-Efficient Simulation
% Load fitted parameters from WhiteNoise_ONOFFalpha_Comparison.m

load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_' stim_wn_id '_Retina_1_MovingNoise_1.mat'], 'OLED');
pix2um = OLED.pixelSize;
% Create spatial filter from fitted parameters
spatial_params = gauss_est(cell_idx, 1:5);
spatial_params(1) = D1_mat/2;
spatial_params(2) = D2_mat/2;
[x, y] = meshgrid(1:D1_mat, 1:D2_mat);
spatial_filter = gaussian2d_one(x, y, spatial_params);
spatial_filter_center = spatial_filter - median(spatial_filter(:)); % Zero mean
spatial_filter_center = spatial_filter_center / max(abs(spatial_filter_center(:))); % Normalize

switch lower(surround_sf_type)
    case 'fixed'
        spatial_params(3:4) = spatial_params(3:4)*4;
    case 'linearfit'
        spatial_params(3:4) = gauss_est(cell_idx, 8:9);
end
spatial_filter = gaussian2d_one(x, y, spatial_params);
spatial_filter_surround = spatial_filter - median(spatial_filter(:));
spatial_filter_surround = spatial_filter_surround / max(abs(spatial_filter_surround(:))); % Normalize
assert(size(spatial_filter_surround, 1) == D2_mat && size(spatial_filter_surround, 2) == D1_mat);
% Create temporal filter from fitted parameters
temporal_params = Gauss_TF_est(cell_idx, :);
temporal_filter = diffgauss_tf(1:D3_mat, temporal_params);
temporal_filter = temporal_filter / sum(abs(temporal_filter)); % Normalize

% spatiotemporal_filter = spatial_filter .* reshape(temporal_filter, 1, 1, []);
STAmat_center = spatial_filter_center .* reshape(temporal_filter, 1, 1, []);
STAmat_center = permute(STAmat_center, [2, 1, 3]);

STAmat_surround = spatial_filter_surround .* reshape(temporal_filter, 1, 1, []);
STAmat_surround = permute(STAmat_surround, [2, 1, 3]);
%%
dr_id = find(dim1_moving_direction == 0);
a2d = @(x) x*pi/180;
[X, Y] = meshgrid(1:D2_mat, 1:D1_mat);
max_time_points = size(Data, 6);
resp = nan(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), max_time_points);
resp_s = nan(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), max_time_points);
for q = 1:length(dim2_contrast)
    for k = 1:length(dim3_bar_width)
        for j = 1:length(dim4_speeds)
            % Initialize stimulus history buffer for temporal filtering
            switch lower(bar_type)
                case 'on'
                    stim_history = (-1 + dim2_contrast(q)*2)*ones(D1_mat, D2_mat, D3_mat);
                case 'off'
                    stim_history = (1 - dim2_contrast(q)*2)*ones(D1_mat, D2_mat, D3_mat);
            end
            
            bw = dim3_bar_width(k)/pix2um;
            diag = 2*(500)/pix2um;
            cAng = dim1_moving_direction(dr_id);
            step = (dim4_speeds(j)/pix2um).*[cos(a2d(cAng+180)) sin(a2d(cAng+180))];
            posi = [cos(a2d(cAng))*0.5*diag sin(a2d(cAng))*0.5*diag]+...
                [cos(a2d(cAng))*0.5*bw sin(a2d(cAng))*0.5*bw];
            num_time = round(min([(diag+bw)*Fz/abs(step(1))+0.5*Fz, max_time_points]));
            
            % Initialize response arrays
            linear_response = nan(1, num_time);
            firing_rate = nan(1, num_time);
            
            for i = 1:num_time
                dstep = i*step/Fz;
                move = posi+dstep;
                
                % Create bar mask for current position
                switch cAng
                    case 0
                        masked = abs((Y-0.5*(D1_mat+1)-move(1))) < 0.5*bw;
                    case 180
                        masked = abs((Y-0.5*(D1_mat+1)-move(1))) < 0.5*bw;
                end

                firing_rate_center(i) = mean(stim_history .* STAmat_center, 'all');
                firing_rate_surround(i) = mean(stim_history .* STAmat_surround, 'all');
                stim_history(:, :, 1) = [];

                % Create current stimulus frame
                switch lower(bar_type)
                    case 'on'
                        canvas = 2*(double(masked)-0.5);
                        canvas(canvas==-1) = -1 + dim2_contrast(q)*2;
                    case 'off'
                        canvas = 2*(-double(masked)+0.5);
                        canvas(canvas==1) = 1 - dim2_contrast(q)*2;
                end
                
                % Update stimulus history (sliding window)
                stim_history(:, :, end+1) = canvas;
                
                if mod(i, 20) == 1 
                    if is_display
                        figure(1);
                        subplot(2, 2, 1);
                        imagesc(canvas', [-1 1]); colorbar;
                        title('Current Stimulus');

                        subplot(2, 2, 2);
                        imagesc(spatial_filter'); colorbar;
                        title('Spatial Filter');

                        subplot(2, 2, 3);
                        plot(1:D3_mat, temporal_filter);
                        title('Temporal Filter');
                        xlabel('Time steps');

                        subplot(2, 2, 4); hold on
                        plot(1:i, firing_rate_center(1:i), 'b');
                        plot(1:i, firing_rate_surround(1:i), 'r');
                        title('Predicted Firing Rate');
                        xlabel('Time steps');
                        ylabel('Firing Rate');

                        drawnow;
                    end
                    clc
                    fprintf('Progress ... %d/%d, %d/%d, %d/%d, %d/%d, %d/%d \n',  ...
                        ii, num_recording, q, length(dim2_contrast), ...
                        k, length(dim3_bar_width), j, length(dim4_speeds), i, num_time);
                end
                
                
            end
            
            % Store results
            resp(q, k, j, 1:length(firing_rate_center)) = firing_rate_center;
            resp_s(q, k, j, 1:length(firing_rate_surround)) = firing_rate_surround;

        end
    end
end


%% Save results with LN model predictions
save_file_name = sprintf('%s_%s_moving_bar_LN_simulated_%d_%s.mat', recording_name,...
    response_name, cell_idx, bar_type);

save(sprintf('./Results/MovingBar/%s', save_file_name), 'dim1_moving_direction', 'dim2_contrast',...
    'dim3_bar_width', 'dim4_speeds','dim5_repeats', 'dim6_time', 'resp', 'resp_s', ...
    'temporal_filter', 'spatial_params', 'temporal_params', 'cell_idx');