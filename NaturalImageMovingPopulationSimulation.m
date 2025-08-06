clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_name = 'b092324';
response_name = 'a082924';
load_recording_name = [recording_name '01'];
Fz = 100;
Speeds = [400 1000 2000 8000];
SurroundType = 1;
bc_thr = -0.6;
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
%%
keyboard;
%% Calulcate line density 
target_density_list = [1 11.69; 5 58.46; 10 116.9; 15 175.3];
num_unit_h = 3;
num_unit_w = 4;
target_area_size = 1000^2; % um^2
pix2um = OLED.pixelSize;
single_unit = sqrt(target_area_size/(num_unit_h*num_unit_w*pix2um^2));
sim_window_h = round(single_unit*num_unit_h); % pixel
sim_window_w = round(single_unit*num_unit_w); % pixel
x_lim = 0.99*[-sim_window_w sim_window_w]/2;
y_lim = 0.99*[-sim_window_h sim_window_h]/2;
target_num_centers = target_density_list(2, 2);
points = create_hexagonal_centers(x_lim, y_lim, target_num_centers, ...
        'NoiseLevel', 0.001, 'RandSeed', 42);
points = round(points);
blurry_length = 50/pix2um; % 0.5*
% Nonlinear parts can be derived from "PBs", "FRs"
%%
dist_m = sqrt((points(:, 1) - points(:, 1)').^2 +(points(:, 2) - points(:, 2)').^2);
dist_m = min(dist_m(eye(size(points, 1))~= 1), [], 'omitnan');
even_split = size(image, 1)/dist_m;
if mod(even_split, 2) == 1
    even_split = even_split -1;
end
line_centers = (size(image, 1)/2-0.5*even_split*dist_m):dist_m:(size(image, 1)/2+0.5*even_split*dist_m);
num_center = length(line_centers);
%% Get model spatial filter
opt_sf = nan(size(opt_sf, 1), size(opt_sf, 2), num_center);
mod_params = optimal_params;
for i = 1:num_center
    mod_params(1) = size(image, 2)/2;
    mod_params(2) = line_centers(i);
    yrange = 200:600;
    switch SurroundType
        case 1
            
        case 2
            mod_params(8:9) = mod_params(3:4)*4;
            mod_params(10) = -mod_params(7)*0.04;
        case 3
            mod_params(8:9) = mod_params(3:4)*4;
            mod_params(10) = -mod_params(7)*0.1;

    end
    copt_sf = gaussian_multi(mod_params, image, num_gauss);
    opt_sf(:, :, i) = copt_sf-median(copt_sf(:));
end    
clear copt_sf
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

opt_STAmat = nan(length(yrange), size(opt_sf, 1), length(opt_tf), num_center);
for i = 1:num_center
    copt_STAmat = opt_sf(:, yrange, i).*reshape(opt_tf, 1, 1, []);
    opt_STAmat(:, :, :, i) = permute(copt_STAmat, [2, 1, 3]);
end
%%
win_h = size(opt_STAmat, 2);
win_w = size(opt_STAmat, 1);
%% Convolute with natural image
img_file_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\MEA\McGillDataset\MatImages';
file_name = 'FileTable.mat';
load(fullfile(img_file_folder, file_name), 'FileTable');
num_image = size(FileTable, 1);
num_speed = length(Speeds);
save_file_name = sprintf('%s_natimg_linesim_srdtype%d_bcthr%0.2G.mat', recording_name, SurroundType, bc_thr);
save_sim_folder = './Results/NatImgSim';
close all
for i = 2:num_image
    file_name = sprintf('ResizeImg%d_%d.mat', FileTable(i, 1), FileTable(i, 2));
    load(fullfile(img_file_folder, file_name), 'img');
    img = squeeze(mean(double(img), 3))/255;
    img = 2*(img-0.5);
    img = img';
    if ~isempty(bc_thr)
        img = BCmodifiedImg(img, bc_thr, blurry_length);
    end
    img_w = size(img, 1);
    img_h = size(img, 2);
    hids = round((0.5*(img_h-win_h)+1):(0.5*(img_h+win_h)));
    simg = img(:, hids);
    img = [flipud(simg(1:size(opt_STAmat, 1), :)); simg; flipud(simg(end-size(opt_STAmat, 1)+1:end, :))];
    mov_dist = size(img, 1)-size(opt_STAmat, 1);
    for j = 1:num_speed
        mov_step = Speeds(j)/(pix2um*Fz);
        if i == 2 && j == 1
            resp = nan(num_image, num_speed, num_center, 2*Fz+round(mov_dist/(Speeds(end)/(pix2um*Fz))));
            cntr = nan(num_image, num_speed, 2*Fz+round(mov_dist/(Speeds(end)/(pix2um*Fz))));
        end
        
        num_step = round(mov_dist/mov_step);
        x = zeros(2*Fz+num_step, 1); 
        for k = 1:num_step
            x(Fz+k, 1) = round(mov_step*k);
        end
        x(Fz+num_step+1:end) = x(Fz+num_step, 1);
        x(:, 2) = x(:, 1) + size(opt_STAmat, 1)-1;
        x = x+1;
        x(x>size(img, 1)) = size(img, 1);
        cmov = repmat(img(x(1, 1):x(1, 2), :), 1, 1, size(opt_STAmat, 3));
        for k = 1:num_step-1
            
            for m = 1:num_center
                % simulation
                % 1 second of exposure to img
                resp(i, j, m, k) = mean(cmov.*squeeze(opt_STAmat(:, :, :, m)), 'all');
                
                clc
                fprintf('Progress ... %d/%d, %d/%d, %d/%d, %d/%d \n',  i, num_image,...
                    j, num_speed, k, num_step-1, m, num_center);
            end
            cmov(:, :, 1) = [];
            canvas = img(x(k+1, 1):x(k+1, 2), :);
            cmov(:, :, end+1) = canvas;
            if mod(k, 10) == 1
                figure(1);
                imagesc(canvas', [-1 1]); colorbar;
            end
            cntr(i, j, k) = std(cmov, [], 'all');
        end
    end
    save(fullfile(save_sim_folder, save_file_name), 'cntr', 'resp', 'bc_thr', 'blurry_length');
    keyboard;
end
%%
foldername = './Results/Params';
% filename = 'b092324_a082924_moving_bar_optimized_6.mat';
% filename = 'b092324_a082924_moving_bar_optimized_6_opt-non_contrast-std_sigmoid0.mat';
filename = 'b092324_a082924_moving_bar_optimized_6_opt-non_contrast-std_sigmoid0_burry1.mat';
load(fullfile(foldername, filename), 'OptW', 'gain_params');
Data = squeeze(resp(2, 1, :, :));
Data(:, isnan(Data(1,:))) = [];
num_center = size(resp, 3);
keepids = std(Data, [], 1) > 0;
Data = Data(:, keepids);
DataY = nan(size(Data));
for i = 1:num_center
    csim = Data(i, :)*1e6;
    % ctr = 1e-7*ones(size(sim));
    % DataY(i, :) = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], sim).*(1./ctr)*OptW(5)+OptW(6).*ctr;
    %         zeros(1, length(sim))], [], 1);
    DataY(i, :) = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*OptW(5)+OptW(6);
            zeros(1, length(csim))], [], 1);
end
figure; 
subplot(1, 2, 1);
imagesc(Data); colorbar
subplot(1, 2, 2);
imagesc(DataY); colorbar
%% Examine the response
fixed = imresize(Data, [size(simg, 2) size(Data, 2)], 'bilinear');
moving = simg';
fixed = (fixed-mean(fixed, 'all'))/std(fixed, [], 'all');
moving = (moving-mean(moving, 'all'))/std(moving, [], 'all');
figure; 
subplot(1, 2, 1);
imagesc(moving); colorbar
subplot(1, 2, 2);
imagesc(fixed); colorbar
%%
[Opt, Met] = imregconfig('Multimodal');
Opt.Epsilon = 1.5e-6;
Opt.MaximumIterations = 200;
Opt.InitialRadius = 6.25e-4;
RefRgt = imregtform(moving, fixed, 'similarity', Opt, Met);
registered = imwarp(moving, RefRgt, 'OutputView', imref2d(size(fixed)));
% Show the images
figure;
subplot(1,3,1);
imagesc(fixed); colorbar
title('Fixed Image');

subplot(1,3,2);
imagesc(moving); colorbar
title('Moving Image');

subplot(1,3,3);
imagesc(registered); colorbar
title(sprintf('Registered Image: %0.3G', val));
%%
% Initial guess for the parameters [scaling, translationX]
initialParams = [1, 0];

% Optimization options
options = optimset('Display', 'iter', 'PlotFcns', @optimplotfval);
CostF = @(p) sum((reshape(LineTransformationF(p, fixed, moving), [], 1) - fixed(:)).^2, 'all', 'omitnan');
% Run the optimizer
clear optimalParams
[optimalParams_x, fval] = fmincon(CostF, [1    0], [], [], [], [],...
                                          [1 -300],...
                                          [1  300]);
CostF = @(p) 1./corr(reshape(LineTransformationF(p, fixed, moving), [], 1) - fixed(:));
[optimalParams_s, fval] = fmincon(CostF, [1    optimalParams_x(2)], [], [], [], [],...
                                       [0.1 optimalParams_x(2)],...
                                       [1  optimalParams_x(2)]);

optimalParams = [optimalParams_s(1) optimalParams_x(2)];
% Display the results
disp('Optimal Parameters:');
disp(['Scaling: ', num2str(optimalParams(1))]);
disp(['Translation X: ', num2str(optimalParams(2))]);

% Apply the optimal transformation  
tformMatrix = [optimalParams(1), 0, 0; 0, 1, 0; optimalParams(2), 0, 1];
% tformMatrix = [1 0 0; 0 1 0; 200 0 1];
tform = affine2d(tformMatrix);
registered = imwarp(moving, tform, 'OutputView', imref2d(size(fixed)));
val = CostF(optimalParams);
% Show the images
figure;
subplot(1,3,1);
imagesc(fixed); colorbar
title('Fixed Image');

subplot(1,3,2);
imagesc(moving); colorbar
title('Moving Image');

subplot(1,3,3);
imagesc(registered); colorbar
title(sprintf('Registered Image: %0.3G', val));
%% Show the images
figure;
subplot(1,3,1);
imshow(fixed);
title('Fixed Image');

subplot(1,3,2);
imshow(moving);
title('Moving Image');

subplot(1,3,3);
imshow(registered);
title('Registered Image');
