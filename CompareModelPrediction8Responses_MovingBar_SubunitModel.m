clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_name = 'b092324';
response_name = 'a082924';
load_recording_name = [recording_name '01'];
Fz = 100;
Speeds = [400 1000 2000 8000];
bc_sf_scaling = 0.2;
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

optimal_params([2 1]) = size(image)/2;
opt_sf = gaussian_multi(optimal_params, image, num_gauss);
opt_sf = opt_sf;
opt_sf = opt_sf-median(opt_sf(:));
spatial_weight = opt_sf;
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

%% Calulcate line density 
target_territory_area = 180 ;
num_unit_h = 3;
num_unit_w = 3;
target_diameter = 300;
target_area_size = 300^2; % um^2
pix2um = OLED.pixelSize;
single_unit = sqrt(target_area_size/(num_unit_h*num_unit_w*pix2um^2));
sim_window_h = round(single_unit*num_unit_h); % pixel
sim_window_w = round(single_unit*num_unit_w); % pixel
x_lim = 0.99*[-sim_window_w sim_window_w]/2;
y_lim = 0.99*[-sim_window_h sim_window_h]/2;
target_num_centers = target_area_size/target_territory_area;
points = create_hexagonal_centers(x_lim, y_lim, target_num_centers, ...
        'NoiseLevel', 0.001, 'RandSeed', 42);
points = round(points);
rmids = sqrt(sum(points.^2, 2))*pix2um > target_diameter*0.5;
points(rmids, :) = [];
%% Get model spatial filter
num_center = size(points, 1);
opt_sf = nan(size(opt_sf, 1), size(opt_sf, 2), num_center);
mod_params = optimal_params;
for i = 1:num_center
    mod_params([2, 1]) = points(i, [2, 1])+size(image)/2;
    mod_params(3:4) = optimal_params(3:4)*bc_sf_scaling;
    
    copt_sf = gaussian_multi(mod_params, image, num_gauss);
    opt_sf(:, :, i) = copt_sf-median(copt_sf(:));
end    
clear copt_sf


%% Get weight matrix from BCs to RGCs
cids = sub2ind(size(image), round(points(:, 2)+size(image, 1)/2), round(points(:, 1)+size(image, 2)/2));
img = zeros(size(image));
img(cids) = 1;
figure; 
subplot(2, 3, 1)
imagesc(smtstdSTA'); colorbar; hold on
title('original STA SF');
subplot(2, 3, 2);
imagesc(squeeze(sum(opt_sf, 3)));colorbar
title('BC array SF');
subplot(2, 3, 3);
plot(mean(squeeze(sum(opt_sf, 3)), 1));
title('BC array SF (slice)');
subplot(2, 3, 4);
imagesc(img); colorbar
title('BC array center');
weights = spatial_weight(cids);
weights = weights/sum(weights);
img = zeros(size(image));
img(cids) = weights;
subplot(2, 3, 5);
imagesc(img); colorbar
title('BC array center weights');
subplot(2, 3, 6);
imagesc(squeeze(sum(opt_sf.*reshape(weights, 1, 1, []), 3))); colorbar
title('BC array SF summation');
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
xrange = 200:400;
yrange = 300:500;
opt_STAmat = nan(length(yrange), length(xrange), length(opt_tf), num_center);
for i = 1:num_center
    copt_STAmat = opt_sf(xrange, yrange, i).*reshape(opt_tf, 1, 1, []);
    opt_STAmat(:, :, :, i) = permute(copt_STAmat, [2, 1, 3]);
end
%% Simulate BCs response to MovingBar and corresponding RGC responses
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
load([stim_data_folder 'Temporal_AlphaRGC_' response_name '_002_Retina_1_MovingBar.mat'], 'MB_IN');

loadFileName = sprintf('%s_moving_bar_processed.mat', response_name);
load(sprintf('./Results/MovingBar/%s', loadFileName), 'dim1_moving_direction', 'dim2_contrast', 'dim3_bar_width',...
    'dim4_speeds', 'dim5_repeats', 'dim6_time', 'Data', 'cData', 'Moving_end_time_ids');
Fz = 100;
%%
cid = size(opt_STAmat, 4);
figure; 
subplot(1, 2, 1);
imagesc(std(squeeze(opt_STAmat(:, :, :, cid)), [], 3)); colorbar
subplot(1, 2, 2);
plot(squeeze(opt_STAmat(111, 133, :, cid))); colorbar
%%
keyboard;
%%
close all
dr_id = 3;
a2d = @(x) x*pi/180;
cmov = -1*ones([size(opt_STAmat, 1) size(opt_STAmat, 2) size(opt_STAmat, 3)]);
num_time = length(x);
t = (0:num_time-1)/Fz;
resp = nan(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), size(Data, 6), num_center);
cntr = nan(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), size(Data, 6));
dim1 = size(cmov, 1);
dim2 = size(cmov, 2);
diag = 500/pix2um;
[X, Y] = meshgrid(1:size(opt_STAmat, 2), 1:size(opt_STAmat, 1));

max_time_points = size(Data, 6);
for q = 1:length(dim2_contrast)
    for k = 1:length(dim3_bar_width)
        for j = 1:length(dim4_speeds)
            cmov = (-1 + dim2_contrast(q)*2)*ones([size(opt_STAmat, 1) size(opt_STAmat, 2) size(opt_STAmat, 3)]);
            bw = dim3_bar_width(k)/pix2um;
            diag = 2*(500)/pix2um;
            cAng = dim1_moving_direction(dr_id);
            step = (dim4_speeds(j)/pix2um).*[cos(a2d(cAng+180)) sin(a2d(cAng+180))];
            posi = [cos(a2d(cAng))*0.5*diag sin(a2d(cAng))*0.5*diag]+...
                [cos(a2d(cAng))*0.5*bw sin(a2d(cAng))*0.5*bw];
            num_time = round(min([(diag+bw)*Fz/step(1)+0.5*Fz, max_time_points]));
            csig = nan(num_time, num_center);
            ccntr = nan(1, num_time);
            for i = 1:num_time
                dstep = i*step/Fz;
                move = posi+dstep;
                switch cAng
                    case 0
                        masked = abs((Y-0.5*(size(opt_STAmat, 1)+1)-move(1))) < 0.5*bw;
                    case 180
                        masked = abs((Y-0.5*(size(opt_STAmat, 1)+1)-move(1))) < 0.5*bw;
                end
                for v = 1:num_center
                    csig(i, v) = mean(cmov.*squeeze(opt_STAmat(:, :, :, v)), 'all');
                end
                ccntr(i) = std(cmov, [], 'all');
                cmov(:, :, 1) = [];
                canvas = 2*(double(masked)-0.5);
                canvas(canvas==-1) = -1 + dim2_contrast(q)*2;
                
                cmov(:, :, end+1) = canvas;
                if mod(i, 10) == 1
                    figure(1);
                    imagesc(canvas', [-1 1]); colorbar;
                end
                clc

                fprintf('Progress ... %d/%d, %d/%d, %d/%d, %d/%d \n',  q, length(dim2_contrast),... 
                    k, length(dim3_bar_width), j, length(dim3_bar_width), i, num_time);
            end
            cntr(q, k, j, 1:size(csig, 1)) = ccntr;
            resp(q, k, j, 1:size(csig, 1), :) = csig;
        end
    end
end
%%
save_file_name = sprintf('%s_%s_moving_bar_simulated_BCModel_sf%0.2G.mat', recording_name,...
        response_name, bc_sf_scaling);
%%
if ~exist(sprintf('./Results/MovingBar/%s', save_file_name), 'file')
    save(sprintf('./Results/MovingBar/%s', save_file_name), 'dim1_moving_direction', 'dim2_contrast',...
        'dim3_bar_width', 'dim4_speeds','dim5_repeats', 'dim6_time', 'resp', 'cntr', 'points', 'pix2um',...
        'weights', 'bc_sf_scaling');
else
    keyboard;
end
%% Get distance dependent exponential relationship between BCs
close all; clc;
lamda = 11.4;
couple_g = -9;
Dists = sqrt((points(:, 1)-points(:, 1)').^2 + (points(:, 2)-points(:, 2)').^2)*pix2um;

dist_decay = exp(-Dists/lamda);
dist_decay(eye(length(dist_decay))==1) = 0;
dist_decay = dist_decay./sum(dist_decay, 2);
dist_decay = -couple_g*dist_decay;
dist_decay(eye(length(dist_decay))==1) = couple_g;



%% Coupling
num_time = size(resp, 4);
resp_cpl_ps = zeros(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), num_time, num_center);
for q = 1:length(dim2_contrast)
    for k = 1:length(dim3_bar_width)
        for j = 1:length(dim4_speeds)
            csig = zeros(num_time, num_center);
            for i = 1:num_time
                t0_sig = squeeze(resp(q, k, j, i, :));
                t0_sig = t0_sig' + t0_sig'*dist_decay';
                csig(i, :) = t0_sig;
            end
            resp_cpl_ps(q, k, j, :, :) = csig;
        end
    end
end

figure; imagesc(squeeze(max(resp_cpl_ps(1, :, 1, :, :), [], [4]))); colorbar
%%
mth_typ = 0;
num_time = size(resp, 4);
resp_cpl = zeros(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), num_time);
for q = 1:length(dim2_contrast)
    for k = 1:length(dim3_bar_width)
        for j = 1:length(dim4_speeds)
            csig = squeeze(resp_cpl_ps(q, k, j, :, :));

            switch mth_typ
                case 0
                case 1
                    csig_r = csig-csig(1, :);
                    csig_r = csig_r.*(max(csig_r, [], 'all')./max(csig_r, [], 1));
                    csig_r = csig_r + csig(1, :);
                    csig = csig_r;
                case 2
                    csig_c = squeeze(resp(q, k, j, :, :));
                    csig_r = csig-csig(1, :);
                    csig_c = csig_c-csig_c(:, 1);
                    csig_r = csig_r.*max(csig_c, [], 1)./max(csig_r, [], 1);
                    csig_r = csig_r + csig(1, :);
                    csig = csig_r;
            end
            resp_cpl(q, k, j, :) = csig*weights;
        end
    end
end

%% normalization
% (1) Normalize to range
% resp_cpl = resp_cpl-min
% resp_cpl = resp_cpl./max(resp_cpl, [], [3 4]);
% (2) Rectified
% resp_cpl = max(resp_cpl, zeros(size(resp_cpl)));
% (3) Renormalize to 0 and 1
% resp_cpl = resp_cpl./max(resp_cpl, [], [3 4]);

%% without coupling
num_time = size(resp, 4);
resp_gwo = zeros(length(dim2_contrast), length(dim3_bar_width), length(dim4_speeds), num_time);
for q = 1:length(dim2_contrast)
    for k = 1:length(dim3_bar_width)
        for j = 1:length(dim4_speeds)
            csig = squeeze(resp(q, k, j, :, :));
            resp_gwo(q, k, j, :) = csig*weights;
        end
    end
end
%%
Linear_blr = load('Z:\Active\Emily\PreyCaptureRGC\Results\MovingBar\b092324_a082924_moving_bar_simulated_6_burry11.4.mat');
Linear_blr = Linear_blr.resp;
Linear_blr = Linear_blr.*max(resp_cpl, [], [2 3 4])./max(Linear_blr, [], [2 3 4]);
Linear = load('Z:\Active\Emily\PreyCaptureRGC\Results\MovingBar\b092324_a082924_moving_bar_simulated_6.mat');
Linear = Linear.resp;
Linear = Linear.*max(resp_cpl, [], [2 3 4])./max(Linear, [], [2 3 4]);
q = 3;
ct = (0:num_time-1)/Fz;
figure;
for k = 1:length(dim3_bar_width)
    for j = 1:length(dim4_speeds)
        subplot(length(dim3_bar_width), length(dim4_speeds),...
            (k-1)*length(dim4_speeds)+j);hold on
        plot(ct, squeeze(resp_cpl(q, k, j, :)), 'Color', [0 0 1]);
        % ylim([0 8e-4])
        % yyaxis right
        plot(ct, squeeze(resp_gwo(q, k, j, :)), 'k');
        
        csig = squeeze(Linear_blr(q, k, j, :));
        plot(ct, csig, 'Color', [1 0 0]);
        csig = squeeze(Linear(q, k, j, :));
        plot(ct, csig, 'Color', 0.5*[0 1 1]);
        xlim([ct(1) 4]);
        xlabel('Time (s)');
        % ylabel('Responses (arbi.)');
        if j == 1
            ylabel(sprintf('Bar width: %d (um)', dim3_bar_width(k)));
            if k == 1
                legend({'Coupling', 'wo coupling', 'image blurred', 'linear'})
            end
        end
        if k == 1
            title(sprintf('Speed: %d (um/s)', dim4_speeds(j)));
        end
        % ylim([-2e-4 8e-4]);
    end
end
sgtitle(sprintf('g: %0.3G bc-sf-scalar: %0.3G', couple_g, bc_sf_scaling))