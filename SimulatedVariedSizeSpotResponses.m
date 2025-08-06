clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_name = 'b092324';
response_name = 'a082924';
load_recording_name = [recording_name '01'];
Fz = 100;
Speeds = [400 1000 2000 8000];
SurroundType = 3;
%% load spot parameters
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
spot_stim_file_name = 'a081024';
load([stim_data_folder 'Temporal_AlphaRGC_' spot_stim_file_name '_001_Retina_1_VariedSizeSpot.mat'], 'SS_IN');
spot_size = SS_IN.diameters;

%% get the linear and nonlinear parts
load([load_data_folder load_recording_name '.mat'], 'masked_STAmat', 'PBs', 'FRs', 'stdSTA');
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_MovingNoise_1.mat'], 'OLED');
pix2um = OLED.pixelSize;
smtstdSTA = medfilt2(stdSTA);
[x, y, dist2d] = peak_distance(smtstdSTA);
minD = 100;
ythr = quantile(y(x>minD), 0.99);
binary_image = smtstdSTA>ythr;
[largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
fprintf('%0.3G', sqrt(sum(largest_segment_mask(:))*4/pi)*pix2um);
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


%% Get model spatial filter
pix2um = OLED.pixelSize;
mod_params = optimal_params;
mod_params(1) = size(image, 2)/2;
mod_params(2) = size(image, 1)/2;
switch SurroundType
    case 1
        yrange = 200:600;
    case 2
        mod_params(8:9) = mod_params(3:4)*4;
        mod_params(10) = -mod_params(7)*0.04;
        yrange = 200:600;
    case 3
        mod_params(8:9) = mod_params(3:4)*4;
        mod_params(10) = -mod_params(7)*0.1;
        yrange = 200:600;

end
opt_sf = gaussian_multi(mod_params, image, num_gauss);
opt_sf = opt_sf-median(opt_sf(:));
figure; 
subplot(1, 2, 1);
imagesc(squeeze(sum(opt_sf, 3)));colorbar
subplot(1, 2, 2);
plot(mean(squeeze(sum(opt_sf, 3)), 1));
%%
keyboard;
%% Get model temporal filter
smtstdSTA = medfilt2(stdSTA);
masked_stdSTA = largest_segment_mask'.*smtstdSTA';
smtSTAmat = medfilt3(masked_STAmat);
smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
tRF = reshape(smtSTAmat, [], size(masked_STAmat, 3))'*masked_stdSTA(:);

[OptimizedParams, scaling] = GaussianTemporalFilter(tRF');
opt_tf = gaussian_temporalfilter(1:length(tRF), OptimizedParams);
opt_tf = opt_tf./sum(abs(opt_tf));

%% Generate 3D model with trim
opt_STAmat = opt_sf(:, yrange).*reshape(opt_tf, 1, 1, []);
opt_STAmat = permute(opt_STAmat, [2, 1, 3]);

%%
diameters = SS_IN.diameters;
diame_ids = 1:length(diameters);
[X, Y] = meshgrid(1:size(opt_STAmat, 2), 1:size(opt_STAmat, 1));
gray_mask = sqrt((X-0.5*(size(opt_STAmat, 2)+1)).^2 +  (Y-0.5*(size(opt_STAmat, 1)+1)).^2) < 0.5*SS_IN.MaskSize/pix2um;
close all
for j = 1:length(diame_ids)
    masked = sqrt((X-0.5*(size(opt_STAmat, 2)+1)).^2 +  (Y-0.5*(size(opt_STAmat, 1)+1)).^2) < 0.5*diameters(j)/pix2um;
    x = zeros(1, Fz*13.5);
    x(51:200) = 1;
    x(201:350) = -1;
    x(351:500) = 1;
    x(501:650) = -1;
    x(651:800) = 1;
    x(801:950) = -1;
    x(951:1100) = 1;
    x(1101:1250) = -1;
    canvas = -1*ones(size(masked));
    canvas(gray_mask == 1) = 0;
    cmov = repmat(canvas, 1, 1, size(opt_STAmat, 3));
    num_time = length(x);
    t = (0:num_time-1)/Fz;
    if j == 1
        resp_ON1 = nan(length(diame_ids), num_time);
    end
    dim1 = size(cmov, 1);
    dim2 = size(cmov, 2);
    for i = 1:num_time
        resp_ON1(j, i) = mean(cmov.*opt_STAmat, 'all');
        cmov(:, :, 1) = [];
        canvas = -1*ones(size(masked));
        canvas(gray_mask == 1) = 0;
        canvas(masked == 1) = x(i);
        cmov(:, :, end+1) = canvas;
        % if mod(i, 10) == 1
        %     figure(1);
        %     imagesc(canvas, [-1 1]); colorbar;
        % end
        clc
        fprintf('Progress ON1... %d/%d, %d/%d \n', j, length(diame_ids), i, num_time);
    end
end
for j = 1:length(diame_ids)
    masked = sqrt((X-0.5*(size(opt_STAmat, 2)+1)).^2 +  (Y-0.5*(size(opt_STAmat, 1)+1)).^2) < 0.5*diameters(j)/pix2um;
    x = zeros(1, Fz*13.5);
    x(51:200) = -1;
    x(201:350) = 1;
    x(351:500) = -1;
    x(501:650) = 1;
    x(651:800) = -1;
    x(801:950) = 1;
    x(951:1100) = -1;
    x(1101:1250) = 1;
    canvas = -1*ones(size(masked));
    canvas(gray_mask == 1) = 0;
    cmov = repmat(canvas, 1, 1, size(opt_STAmat, 3));
    num_time = length(x);
    t = (0:num_time-1)/Fz;
    if j == 1
        resp_OFF1 = nan(length(diame_ids), num_time);
    end
    dim1 = size(cmov, 1);
    dim2 = size(cmov, 2);
    for i = 1:num_time
        resp_OFF1(j, i) = mean(cmov.*opt_STAmat, 'all');
        cmov(:, :, 1) = [];
        canvas = -1*ones(size(masked));
        canvas(gray_mask == 1) = 0;
        canvas(masked == 1) = x(i);
        cmov(:, :, end+1) = canvas;
        % if mod(i, 10) == 1
        %     figure(1);
        %     imagesc(canvas, [-1 1]); colorbar;
        % end
        clc
        fprintf('Progress OFF1... %d/%d, %d/%d \n', j, length(diame_ids), i, num_time);
    end
end
%%
save_file_name = sprintf('%s_stationaryspot_simulated_srdtype%d.mat', recording_name, ...
    SurroundType);
%%
save(sprintf('./Results/VariedSizeSpotSim/%s', save_file_name), 'diame_ids', 'x',...
    'resp_ON1', 'resp_OFF1');

%%
colors = parula (4);
figure; 
for i = 1:4
    subplot(1, 2, 1); hold on
    plot(squeeze(Data(7, 2, i, :)), 'Color', colors(i, :));
    subplot(1, 2, 2); hold on
    plot(squeeze(Data(7, 2, i+4, :)), 'Color', colors(i, :));
end
%%
cuts = [301 601 901];
time_span = 2*Fz;
sum_resp = nan(size(resp_ON1, 1), 1);
for i = 1:size(resp_ON1, 1)
    ctrace = [];
    for j = 1:3
        ctrace = [ctrace; resp_ON1(i, cuts(j):cuts(j)+time_span-1)];
    end
    sum_resp(i) = mean(ctrace(:, 51:80) - ctrace(:, 1:30), 'all');
end
figure; plot(diameters, sum_resp);
ylim([0 1.2*max(sum_resp)]);
xlim([0 diameters(end)]);
box off
ylabel('Resp (arbi.)');
xlabel('Diameters (um)');
title(sprintf('%s Surround type: %d', recording_name, SurroundType))