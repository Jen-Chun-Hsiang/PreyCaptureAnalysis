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
blurry_length = 50/pix2um; % 0.5*

%%
num_gauss = 1;
image = stdSTA';
initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];
initial_params = repmat(initial_params, num_gauss, 1);
initial_params(:, 1:2) = initial_params(:, 1:2) + 10*rand(num_gauss, 2);
initial_params = initial_params';
objective_function = @(params) 1-corr(image(:), reshape(gaussian_multi(params, image, num_gauss), [], 1));
options.MaxFunEvals = 600*length(initial_params(:));
[SF_params,fval,exitflag,output]= fminsearch(objective_function, initial_params);

opt_sf = gaussian_multi(SF_params, image, num_gauss);
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
%% Get model temporal filter
smtstdSTA = medfilt2(stdSTA);
masked_stdSTA = largest_segment_mask'.*smtstdSTA';
smtSTAmat = medfilt3(masked_STAmat);
smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
tRF = reshape(smtSTAmat, [], size(masked_STAmat, 3))'*masked_stdSTA(:);

[OptimizedParams, scaling] = GaussianTemporalFilter(tRF');
%%
TF_params = gaussian_temporalfilter(1:length(tRF), OptimizedParams);

mov_params = OptimizedParams;
mov_tf = gaussian_temporalfilter(1:length(tRF), mov_params);
opt_tf = gaussian_temporalfilter(1:length(tRF), OptimizedParams);
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
% keyboard;
%%
save_file_name = sprintf('%s_LinearNL_params.mat', recording_name);
save(sprintf('./Results/MovingBar/%s', save_file_name), 'SF_params', 'TF_params', 'tRF', 'smtstdSTA');