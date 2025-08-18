%% Modified LN Model Implementation for Memory-Efficient Simulation
% Load fitted parameters from WhiteNoise_ONOFFalpha_Comparison.m
process_version = 'GaussianFitting_processed_081425_1.mat';
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
processedFile = fullfile(folder_name, process_version);

if exist(processedFile, 'file')
    load(processedFile, 'gauss_est', 'Gauss_TF_est', 'NL_curves', 'NL_params');
    fprintf('Loaded fitted parameters for LN model\n');
else
    error('Fitted parameters not found. Run WhiteNoise_ONOFFalpha_Comparison.m first.');
end

% Select cell for simulation (you can change this based on your current cell)
cell_idx = 1; % Modify this to match your current cell being processed

% Create spatial filter from fitted parameters
spatial_params = gauss_est(cell_idx, :);

% modify the center-surround ratio based on the measurement
spatial_filter = gaussian_multi(spatial_params, zeros(size(opt_STAmat, 1), size(opt_STAmat, 2)), 1);
spatial_filter = spatial_filter - median(spatial_filter(:)); % Zero mean

% Create temporal filter from fitted parameters
temporal_params = Gauss_TF_est(cell_idx, :);
temporal_length = size(opt_STAmat, 3);
temporal_filter = gaussian_temporalfilter(1:temporal_length, temporal_params);
temporal_filter = temporal_filter / sum(abs(temporal_filter)); % Normalize

% Create nonlinear function
if exist('NL_params', 'var') && size(NL_params, 1) >= cell_idx
    nl_params = NL_params(cell_idx, :);
    nonlinear_func = @(x) nl_params(1) * normcdf(x, nl_params(2), nl_params(3)) + nl_params(4);
elseif exist('NL_curves', 'var') && size(NL_curves, 1) >= cell_idx
    x_sample_range = linspace(-3, 3, size(NL_curves, 2));
    nl_curve = NL_curves(cell_idx, :);
    nonlinear_func = @(x) interp1(x_sample_range, nl_curve, x, 'linear', 'extrap');
else
    % Fallback: simple rectification
    nonlinear_func = @(x) max(0, x);
end

spatiotemporal_filter = spatial_filter .* reshape(temporal_filter, 1, 1, []);

%%
max_time_points = size(Data, 6);
for q = 1:length(dim2_contrast)
    for k = 1:length(dim3_bar_width)
        for j = 1:length(dim4_speeds)
            % Initialize stimulus history buffer for temporal filtering
            switch lower(bar_type)
                case 'on'
                    stim_history = (-1 + dim2_contrast(q)*2)*ones(size(opt_STAmat, 1), size(opt_STAmat, 2), temporal_length);
                case 'off'
                    stim_history = (1 - dim2_contrast(q)*2)*ones(size(opt_STAmat, 1), size(opt_STAmat, 2), temporal_length);
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
                        masked = abs((Y-0.5*(size(opt_STAmat, 1)+1)-move(1))) < 0.5*bw;
                    case 180
                        masked = abs((Y-0.5*(size(opt_STAmat, 1)+1)-move(1))) < 0.5*bw;
                end
                
                % Create current stimulus frame
                switch lower(bar_type)
                    case 'on'
                        canvas = 2*(double(masked)-0.5);
                        canvas(canvas==-1) = -1 + dim2_contrast(q)*2;
                        if is_blurry
                            canvas = imgaussfilt(canvas, blurry_length);
                            if j == 1 && i == 1
                                maxv = max(canvas(:));
                                minv = min(canvas(:));
                            end
                            extend_ratio = (1-minv)/(maxv-minv);
                            canvas = canvas*extend_ratio+minv*(1-extend_ratio);
                        end
                    case 'off'
                        canvas = 2*(-double(masked)+0.5);
                        canvas(canvas==1) = 1 - dim2_contrast(q)*2;
                        if is_blurry
                            canvas = imgaussfilt(canvas, blurry_length);
                            if j == 1 && i == 1
                                maxv = max(canvas(:));
                                minv = min(canvas(:));
                            end
                            extend_ratio = (-1-maxv)/(minv-maxv);
                            canvas = canvas*extend_ratio+maxv*(1-extend_ratio);
                        end
                end
                
                % Update stimulus history (sliding window)
                stim_history(:, :, 1) = [];
                stim_history(:, :, end+1) = canvas;
                
                % Compute linear response using LN model
                % 1. Spatial filtering for each time point in history
                linear_response(i) = mean(stim_history .* spatiotemporal_filter, 'all');
                
                % 3. Apply nonlinearity
                firing_rate(i) = nonlinear_func(linear_response(i));
                
                if mod(i, 10) == 1 && is_display
                    figure(1);
                    subplot(2, 2, 1);
                    imagesc(canvas', [-1 1]); colorbar;
                    title('Current Stimulus');
                    
                    subplot(2, 2, 2);
                    imagesc(spatial_filter'); colorbar;
                    title('Spatial Filter');
                    
                    subplot(2, 2, 3);
                    plot(1:temporal_length, temporal_filter);
                    title('Temporal Filter');
                    xlabel('Time steps');
                    
                    subplot(2, 2, 4);
                    plot(1:i, firing_rate(1:i));
                    title('Predicted Firing Rate');
                    xlabel('Time steps');
                    ylabel('Firing Rate');
                    
                    drawnow;
                end
                
                fprintf('Progress ... %d/%d, %d/%d, %d/%d, %d/%d, %d/%d, %d/%d \n',  ...
                    ii, num_recording, jj, 2, q, length(dim2_contrast), ...
                    k, length(dim3_bar_width), j, length(dim4_speeds), i, num_time);
            end
            
            % Store results
            resp(q, k, j, 1:length(firing_rate)) = firing_rate;
            resp_linear(q, k, j, 1:length(linear_response)) = linear_response;
            
            % Optional: store some diagnostics
            if q == 1 && k == 1 && j == 1
                % Save example temporal trace for debugging
                example_spatial_trace = spatial_responses;
                example_linear_trace = linear_response;
                example_nl_trace = firing_rate;
            end
        end
    end
end


%% Save results with LN model predictions
if is_blurry
    save_file_name = sprintf('%s_%s_moving_bar_LN_simulated_%d_burry%0.3G_%s.mat', recording_name,...
        response_name, cell_idx, blurry_length, bar_type);
else
    save_file_name = sprintf('%s_%s_moving_bar_LN_simulated_%d_%s.mat', recording_name,...
        response_name, cell_idx, bar_type);
end

save(sprintf('./Results/MovingBar/%s', save_file_name), 'dim1_moving_direction', 'dim2_contrast',...
    'dim3_bar_width', 'dim4_speeds','dim5_repeats', 'dim6_time', 'resp', 'resp_linear', ...
    'spatial_filter', 'temporal_filter', 'spatial_params', 'temporal_params', 'cell_idx', ...
    'example_spatial_trace', 'example_linear_trace', 'example_nl_trace');