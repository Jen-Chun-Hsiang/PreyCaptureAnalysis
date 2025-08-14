%% get the linear and nonlinear parts
load([load_data_folder load_recording_name '.mat'], 'masked_STAmat', 'PBs', 'FRs', 'stdSTA');
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_' stim_wn_id '_Retina_1_MovingNoise_1.mat'], 'OLED');
pix2um = OLED.pixelSize;
smtstdSTA = medfilt2(stdSTA);
[x, y, dist2d] = peak_distance(smtstdSTA);
minD = 100;
ythr = quantile(y(x>minD), 0.99);
binary_image = smtstdSTA>ythr;
[largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
fprintf('%0.3G', sqrt(sum(largest_segment_mask(:))*4/pi)*pix2um);
%%
if is_blurry
    blurry_length = 50/pix2um; % 0.5*
end
%%
num_gauss = 1;
image = stdSTA';
initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];
initial_params = repmat(initial_params, num_gauss, 1);
initial_params(:, 1:2) = initial_params(:, 1:2) + 10*rand(num_gauss, 2);
initial_params = initial_params';
objective_function = @(params) 1-corr(image(:), reshape(gaussian_multi(params, image, num_gauss), [], 1));
options.MaxFunEvals = 600*length(initial_params(:));
[optimal_params,fval,exitflag,output]= fminsearch(objective_function, initial_params);

opt_sf = gaussian_multi(optimal_params, image, num_gauss);
opt_sf = opt_sf;
opt_sf = opt_sf-median(opt_sf(:));
figure; 
subplot(1, 2, 1)
imagesc(smtstdSTA'); colorbar; hold on
plot([400 400], [0 600], '--k');
plot([0 800], [300 300], '--k');
subplot(1, 2, 2)
imagesc(opt_sf); colorbar; hold on
plot([400 400], [0 600], '--k');
plot([0 800], [300 300], '--k');
% Spatial temporal filter is in "masked_STAmat"

%%
% keyboard;
% Nonlinear parts can be derived from "PBs", "FRs"
%% Get model spatial filter
optimal_params(1) = size(image, 2)/2;
optimal_params(2) = size(image, 1)/2;
opt_sf = gaussian_multi(optimal_params, image, num_gauss);
opt_sf = opt_sf-median(opt_sf(:));



%% 

mov_sf_params = optimal_params;
switch implement_case_id
    case {1, 2, 3}
        mov_sf = gaussian_multi(mov_sf_params, image, num_gauss);
        mov_sf = mov_sf-median(mov_sf(:));
    case 4
        mov_sf_params(3:4) = mov_sf_params(3:4)*0.7;
        mov_sf = gaussian_multi(mov_sf_params, image, num_gauss);
        mov_sf = mov_sf-median(mov_sf(:));
        mov_sf = mov_sf-2*mean(mov_sf(:)); % from 50
    case 5
        mov_sf_params(3:4) = mov_sf_params(3:4)*0.5;
        mov_sf = gaussian_multi(mov_sf_params, image, num_gauss);
        mov_sf = mov_sf-median(mov_sf(:));
        mov_sf = mov_sf-2*mean(mov_sf(:)); % from 50
    case 6
        mov_sf_params(3:4) = mov_sf_params(3:4)*0.5;
        mov_sf = gaussian_multi(mov_sf_params, image, num_gauss);
        mov_sf = mov_sf-median(mov_sf(:));
    case 7
        mov_sf_params(3:4) = mov_sf_params(3:4)*0.5;
        mov_sf = gaussian_multi(mov_sf_params, image, num_gauss);
        mov_sf = mov_sf-median(mov_sf(:));
        optimal_params(3:4) = optimal_params(3:4)*1.5;
        opt_sf = gaussian_multi(optimal_params, image, num_gauss);
        opt_sf = opt_sf-median(opt_sf(:));

end
figure; 
subplot(1, 2, 1)
imagesc(opt_sf); colorbar
title('Excitation spatial filter');
subplot(1, 2, 2);
imagesc(mov_sf); colorbar
title('Inhibition spatial filter');
%%
% keyboard;
%% Get model temporal filter
smtstdSTA = medfilt2(stdSTA);
masked_stdSTA = largest_segment_mask'.*smtstdSTA';
smtSTAmat = medfilt3(masked_STAmat);
smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
tRF = reshape(smtSTAmat, [], size(masked_STAmat, 3))'*masked_stdSTA(:);

[OptimizedParams, scaling] = GaussianTemporalFilter(tRF');
%%
opt_tf = gaussian_temporalfilter(1:length(tRF), OptimizedParams);

switch implement_case_id
    case 1
        delay_fac = 7;
        mov_params = OptimizedParams;
        mov_params(3:4) = mov_params(3:4)-delay_fac;
    case 2
        delay_fac = 4;
        mov_params = OptimizedParams;
        mov_params(3:4) = mov_params(3:4)-delay_fac;
        mov_params(6) = 0.1;
    case 3
        delay_fac = 4;
        mov_params = OptimizedParams;
        mov_params(3:4) = mov_params(3:4)-delay_fac;
        mov_params(6) = 0.1;
        mov_params(1:2) = mov_params(1:2)*1.5;
    case {4, 5, 6, 7}
        mov_params = OptimizedParams;
        
end

mov_tf = gaussian_temporalfilter(1:length(tRF), mov_params);

opt_tf = opt_tf./sum(abs(opt_tf));
mov_tf = mov_tf./sum(abs(mov_tf));
tRF = tRF./sum(abs(tRF));

figure; hold on
t = -0.5:1/Fz:0;
t(end) = [];
plot(t, tRF, 'k');
plot(t, opt_tf, 'b');
plot(t, mov_tf, 'm');
xlabel('Time (s)');
ylabel('Effective contrast');
legend({'TF', 'exc TF', 'inh TF'});
%% Generate Combine STA
opt_STAmat = opt_sf.*reshape(opt_tf, 1, 1, []);
opt_STAmat = permute(opt_STAmat, [2, 1, 3]);
masked_STAmat = opt_STAmat;
mov_STAmat = mov_sf.*reshape(mov_tf, 1, 1, []);
mov_STAmat = permute(mov_STAmat, [2, 1, 3]);

%%
% keyboard;
%%
[nl_fuc, divider]  = getNonlinearFunc(PBs, FRs);
%% Simulation (2) moving bars
BarWidth = 200;
Speed = 2000;
cAng = 0;
a2d = @(x) x*pi/180;
sw = OLED.width;
sh = OLED.height;
diag = sqrt(sw^2+sh^2);
bw = BarWidth/pix2um;
posi = [0 0];
step = Speed.*[cos(a2d(cAng+180)) sin(a2d(cAng+180))];
[X, Y] = meshgrid(1:size(opt_STAmat, 2), 1:size(opt_STAmat, 1));
switch cAng
    case 0
        masked = abs((Y-0.5*(size(opt_STAmat, 1)+1)-posi(1))) < 0.5*bw;
end

%%
figure;
subplot(1, 3, 1)
imagesc(masked');
subplot(1, 3, 2)
imagesc(opt_STAmat(:, :, 45)');
subplot(1, 3, 3)
imagesc(mov_STAmat(:, :, 45)');

%%
% posi = [cos(a2d(cAng))*0.5*diag sin(a2d(cAng))*0.5*diag]+...
%     [cos(a2d(cAng))*0.5*bw sin(a2d(cAng))*0.5*bw];
% x = ones(1, Fz*3.5);
% cmov = -1*ones(size(opt_STAmat));
% num_time = length(x);
% t = (0:num_time-1)/Fz;
% resp = nan(1, num_time);
% dim1 = size(cmov, 1);
% dim2 = size(cmov, 2);
% for i = 1:num_time
%     dstep = i*step/Fz;
%     move = posi+dstep;
%     switch cAng
%         case 0
%             masked = abs((Y-0.5*(size(opt_STAmat, 1)+1)-move(1))) < 0.5*bw;
%     end
%     resp(i) = mean(cmov.*opt_STAmat, 'all') - mean(cmov.*mov_STAmat, 'all');
%     cmov(:, :, 1) = [];
%     cmov(:, :, end+1) = 2*(double(masked)-0.5);
%     clc
%     fprintf('Progress ... %d/%d \n', i, num_time);
% end
%%
% figure;
% subplot(1, 2, 1)
% plot(t, resp, 'k');
% subplot(1, 2, 2)
% % plot(t, nl_fuc(resp*0.1/divider), 'k');
% plot(t, resp, 'k');

%% Simulation (1) varied size spots
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
load([stim_data_folder 'Temporal_AlphaRGC_' response_name '_' stim_mb_id '_Retina_1_MovingBar.mat'], 'MB_IN');

loadFileName = sprintf('%s_moving_bar_processed.mat', response_name);
load(sprintf('./Results/MovingBar/%s', loadFileName), 'dim1_moving_direction', 'dim2_contrast', 'dim3_bar_width',...
    'dim4_speeds', 'dim5_repeats', 'dim6_time', 'Data', 'cData', 'Moving_end_time_ids');
Fz = 100;
%%
% keyboard;
%%
dr_id = 2;
ct_id = 1;
bw_id = 3;
sp_id = 1;
figure;
% ct = (0:round(mean(Moving_end_time_ids(dr_id, bw_id, sp_id, :), 4, 'omitnan'))-1)/Fz
ct = (0:size(Data, 6)-1)/Fz;
plot(ct, squeeze(mean(Data(dr_id, ct_id, bw_id, sp_id, :, :), 5, 'omitnan')));
title(sprintf('Contrast:%0.2G %d (um) %d (um/s)', dim2_contrast(ct_id), dim3_bar_width(bw_id), dim4_speeds(sp_id)));
%%
dr_id = 3;
%%
close all
num_time = length(x);
t = (0:num_time-1)/Fz;
resp = nan(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), size(Data, 6));
resp_s = nan(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), size(Data, 6));
cntr = nan(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), size(Data, 6));
dim1 = size(opt_STAmat, 1);
dim2 = size(opt_STAmat, 2);
diag = 500/pix2um;

max_time_points = size(Data, 6);
[Y, X] = meshgrid(1:size(opt_STAmat, 1), 1:size(opt_STAmat, 2)); % Precompute meshgrid

for q = 1:length(dim2_contrast)
    for k = 1:length(dim3_bar_width)
        bw = dim3_bar_width(k)/pix2um;
        for j = 1:length(dim4_speeds)
            % Precompute constants outside time loop
            diag = 2*(500)/pix2um;
            cAng = dim1_moving_direction(dr_id);
            step = (dim4_speeds(j)/pix2um).*[cos(a2d(cAng+180)) sin(a2d(cAng+180))];
            posi = [cos(a2d(cAng))*0.5*diag sin(a2d(cAng))*0.5*diag] + ...
                   [cos(a2d(cAng))*0.5*bw sin(a2d(cAng))*0.5*bw];
            num_time = round(min([(diag+bw)*Fz/step(1)+0.5*Fz, max_time_points]));
            csig = nan(1, num_time);
            csig_s = nan(1, num_time);
            ccntr = nan(1, num_time);

            % Initialize cmov only once per (q,k,j)
            switch lower(bar_type)
                case 'on'
                    cmov = (-1 + dim2_contrast(q)*2)*ones([size(opt_STAmat), num_time+1]);
                case 'off'
                    cmov = (1 - dim2_contrast(q)*2)*ones([size(opt_STAmat), num_time+1]);
            end

            % Precompute mask shifts for all time points if possible
            mask_shifts = zeros(size(Y,1), size(Y,2), num_time);
            for i = 1:num_time
                dstep = i*step/Fz;
                move = posi + dstep;
                switch cAng
                    case {0, 180}
                        mask_shifts(:,:,i) = abs((Y-0.5*(size(opt_STAmat, 1)+1)-move(1))) < 0.5*bw;
                end
            end

            for i = 1:num_time
                masked = mask_shifts(:,:,i);

                % Only use the last cmov slice for calculation
                cmov_slice = cmov(:,:,i);

                csig(i) = mean(cmov_slice.*opt_STAmat, 'all');
                csig_s(i) = mean(cmov_slice.*mov_STAmat, 'all');
                ccntr(i) = std(cmov_slice, [], 'all');

                % Generate canvas
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

                % Store canvas in cmov
                cmov(:,:,i+1) = canvas;

                % Only show figure every 50 steps or as needed
                if mod(i, 50) == 1
                    % Uncomment for debugging
                    % figure(1); imagesc(canvas', [-1 1]); colorbar;
                end
            end

            % Assign results
            cntr(q, k, j, 1:num_time) = ccntr;
            resp(q, k, j, 1:num_time) = csig;
            resp_s(q, k, j, 1:num_time) = csig_s;

            % Print progress every outer loop
            fprintf('Progress ... q=%d/%d, k=%d/%d, j=%d/%d\n', ...
                q, length(dim2_contrast), k, length(dim3_bar_width), j, length(dim4_speeds));
        end
    end
end

%% save
if is_blurry
    save_file_name = sprintf('%s_%s_moving_bar_simulated_%d_burry%0.3G_%s.mat', recording_name,...
        response_name, implement_case_id, blurry_length, bar_type);
else
    save_file_name = sprintf('%s_%s_moving_bar_simulated_%d_%s.mat', recording_name,...
        response_name, implement_case_id, bar_type);
end
%%
save(sprintf('./Results/MovingBar/%s', save_file_name), 'dim1_moving_direction', 'dim2_contrast',...
    'dim3_bar_width', 'dim4_speeds','dim5_repeats', 'dim6_time', 'resp', 'resp_s', 'cntr');
