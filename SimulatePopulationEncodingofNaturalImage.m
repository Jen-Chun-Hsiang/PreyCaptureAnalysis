clear all; close all; clc;
%% Want to simulate how density affect image recognition
% [Analysis notes]
% (1) Use PCA to find the low dimension that present the highest variation
% across sample naturalistic images
% (2) Use the projection in step 1 to see how it can distinguish across
% distance moving of the same images
% (3) distinguish score is (distance between image and moved version) /
% (average distance across difference images)
% Hypothesis is that for lower density, it is harder to separate the
% distance of slightly move images
% Plot distinguish score vs moving distance (um)

% [Simulation notes]
% (1) Use the experimental data for the NL model, cut the center of SF and
% move to the simulated location
% (2) Randomly sample from experimental data, one data one time, so repeat
% measurement will be the data paired we have

%%
% Density of simulated ONs-RGC
% 1 (count): 11.69 /mm2 (density), 5: 58.46 /mm2, 10: 116.9/mm2,
% 15:175.3/mm2

%% parameters
sim_data_time = '08222402'; % mark the simulation
exp_cell_name = 'a081024'; % c081224file list to get NLmodel

random_seed = 42;
num_mov_image = 40;
num_mov_sample = 20;
move_distance = 10; % 99
is_display = 0;

% less changed variables
Fz = 100;
target_area_size = 1000^2; % um^2
target_density_list = [1 11.69; 5 58.46; 10 116.9; 15 175.3]; % 30 350.7];
num_unit_h = 3;
num_unit_w = 4;

load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';

%% load data
recording_name = exp_cell_name;
load_recording_name = [recording_name '01'];
load([load_data_folder load_recording_name '.mat'], 'masked_STAmat', 'PBs', 'FRs', 'stdSTA', 'tRF');
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_MovingNoise_1.mat'], 'OLED');

pix2um = OLED.pixelSize;
smtstdSTA = medfilt2(stdSTA);
[x, y, dist2d] = peak_distance(smtstdSTA);
minD = 100;
ythr = quantile(y(x>minD), 0.99);
binary_image = smtstdSTA>ythr;
[largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
mask_r_pix = sqrt(sum(largest_segment_mask(:))/pi);
mask_r = mask_r_pix*pix2um;
fprintf('%0.3G (diameters)', mask_r*2);

[~, tf_peak_id] = max(tRF);

%%
num_gauss = 1;
image = smtstdSTA';
initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];
initial_params = repmat(initial_params, num_gauss, 1);
initial_params(:, 1:2) = initial_params(:, 1:2) + 10*rand(num_gauss, 2);
initial_params = initial_params';
objective_function = @(params) 1-corr(image(:), reshape(gaussian_multi(params, image, num_gauss), [], 1));
options.MaxFunEvals = 600*length(initial_params(:));
[optimal_params,fval,exitflag,output]= fminsearch(objective_function, initial_params);

figure;
subplot(1, 2, 1);
imagesc(smtstdSTA); colorbar; hold on;
plot(optimal_params(2), optimal_params(1), 'rx');
subplot(1, 2, 2);
imagesc(largest_segment_mask); colorbar;
%%
modified_params = optimal_params;
surround_fac = 3;
modified_params(10) =  -1;
modified_params(8) = optimal_params(3)*surround_fac;
modified_params(9) = optimal_params(4)*surround_fac;
gaussian_model_optimal = gaussian_multi(optimal_params, image, num_gauss);
gaussian_model_modified = gaussian_multi(modified_params, image, num_gauss);
figure; 
subplot(1, 2, 1)
imagesc(gaussian_model_optimal); colorbar
subplot(1, 2, 2)
imagesc(gaussian_model_modified); colorbar

%%
clear Data
num_density = length(target_density_list);
single_unit = sqrt(target_area_size/(num_unit_h*num_unit_w*pix2um^2));
sim_window_h = round(single_unit*num_unit_h); % pixel
sim_window_w = round(single_unit*num_unit_w); % pixel
x_lim = 0.99*[-sim_window_w sim_window_w]/2;
y_lim = 0.99*[-sim_window_h sim_window_h]/2;

win_h = ceil(sim_window_h+2*mask_r_pix);
win_w = ceil(sim_window_w+2*mask_r_pix);

%generate an array of cut receptive field
cube_STA = extract_STA_center(masked_STAmat, optimal_params(1:2), mask_r_pix);
cube_x = size(cube_STA, 1);
cube_y = size(cube_STA, 2);


for j = 4:num_density
    target_density_id = j;
    %% estimation of sizes
    target_num_centers = target_density_list(target_density_id, 2);
    points = create_hexagonal_centers(x_lim, y_lim, target_num_centers, ...
        'NoiseLevel', 0.3, 'RandSeed', 42);
    points = round(points);

    Data{j, 3} = points;
    %%
    is_density_map = 0;
    if is_density_map
        density_map = zeros(400, 300);
        [X, Y] = meshgrid(1:400, 1:300);
        for i = 1:size(points, 1)
            map_ids = sqrt((X'-200+1-points(i, 1)).^2 + (Y'-150+1-points(i, 2)).^2) < mask_r_pix;
            density_map = density_map+double(map_ids);
        end
        figure;
        subplot(1, 2, 1), hold on
        scatter(points(:, 1), points(:, 2), 5, 'k', 'filled');
        viscircles(points, mask_r_pix, 'Color','m');
        subplot(1, 2, 2);
        imagesc(density_map'); colorbar
        sgtitle(sprintf('Density: %d, Cell count: %d', target_density_list(target_density_id, 1), size(points, 1)));
        keyboard;
    end
    %% Going through all the images and movement of images
    close all
    num_sim_cell = size(points, 1);
    rf_maps = nan(win_w, win_h, num_sim_cell);
    for i = 1:num_sim_cell
        cmap = zeros(win_w, win_h, size(cube_STA, 3));
        hids = (points(i, 2)+round(0.5*win_h-mask_r_pix)+1):(points(i, 2)+round(0.5*win_h+mask_r_pix));
        wids = (points(i, 1)+round(0.5*win_w-mask_r_pix)+1):(points(i, 1)+round(0.5*win_w+mask_r_pix));
        assert(length(hids) == cube_x);
        cmap(wids, hids, :) = cube_STA;

        % figure(1);
        % subplot(1, 2, 1)
        % imagesc(std(cmap, [], 3));
        % subplot(1, 2, 2)
        % imagesc(squeeze(cmap(:, :, tf_peak_id)));
        % pause(1)
        rf_maps(:, :, i) = cmap(:, :, tf_peak_id);
    end
    rf_maps = reshape(rf_maps, [], num_sim_cell);
    rf_maps = rf_maps./sum(rf_maps);
    %%
    img_file_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\MEA\McGillDataset\MatImages';
    file_name = 'FileTable.mat';
    load(fullfile(img_file_folder, file_name), 'FileTable');
    num_image = size(FileTable, 1);
    rng(random_seed)
    moving_test_image_id = randperm(num_image, num_mov_image);
    moving_x = randsample(-move_distance:move_distance, num_mov_sample);
    moving_x = sort(moving_x);
    moving_y = randsampe(-move_distance:move_distance, num_mov_sample);
    moving_y = sort(moving_y);

    cData = nan(num_image, num_sim_cell);
    % every natural image
    for i = 1:num_image
        file_name = sprintf('ResizeImg%d_%d.mat', FileTable(i, 1), FileTable(i, 2));
        load(fullfile(img_file_folder, file_name), 'img');
        img = squeeze(mean(double(img), 3))/255;
        img = 2*(img-0.5);
        img = img';
        img_w = size(img, 1);
        img_h = size(img, 2);
        hids = round((0.5*(img_h-win_h)+1):(0.5*(img_h+win_h)));
        wids = round((0.5*(img_w-win_w)+1):(0.5*(img_w+win_w)));
        img = img(wids, hids);


        % figure(1);
        % imagesc(img, [0 1]);
        cimg = img(:)'*rf_maps;
        cData(i, :) = cimg;
        % cimg = (cimg-min(cimg))/range(cimg);
        close all
        if is_display
            figure(1);
            subplot(1, 2, 1);
            imagesc(img');
            colormap(gray);
            subplot(1, 2, 2); hold on
            
            for k = 1:num_sim_cell
                scatter(points(k, 1), -points(k, 2), 25, cimg(k)*ones(1, 3));
            end
            xlim([-0.5*win_w 0.5*win_w])
            ylim([-0.5*win_h 0.5*win_h])
        end
        dot_points = round(points+size(img)*0.5);
        dot_ids = sub2ind(size(img), dot_points(:, 1), dot_points(:, 2));
        dot_value = img(dot_ids);

        

        rcimg = reconstructImagebyPoints(img, points(:, 1), points(:, 2), dot_value);
        rimg = reconstructImagebyPoints(img, points(:, 1), points(:, 2), cimg);

        img = (img-min(img(:)))/range(img(:));
        rcimg = (rcimg-min(rcimg(:)))/range(rcimg(:));
        rimg = (rimg-min(rimg(:)))/range(rimg(:));
        figure(2); 
        subplot(1, 3, 1);
        imagesc(img'); colorbar;
        subplot(1, 3, 2);
        imagesc(rimg); colorbar;
        subplot(1, 3, 3);
        imagesc(rcimg); colorbar;
        keyboard;
        
        % not yet finished
        [spatail_freq, freq_power] = getSpatialFreqPower(img);
        [spatail_freq_r, freq_power_r] = getSpatialFreqPower(rimg);
        [spatail_freq_rc, freq_power_rc] = getSpatialFreqPower(rcimg);
        if i == 1
            col_freq_power = zeros(num_image, length(freq_power));
            col_freq_power_r= zeros(num_image, length(freq_power_r));
            col_freq_power_rc = zeros(num_image, length(freq_power_rc));
        end
        col_freq_power(i, :) = freq_power;
        col_freq_power_r(i, :) = freq_power_r;
        col_freq_power_rc(i, :) = freq_power_rc;

        if is_display
            figure(3); hold on
            plot(spatail_freq, log10(freq_power), 'k');
            plot(spatail_freq_r, log10(freq_power_r), 'r');
            plot(spatail_freq_rc, log10(freq_power_rc), 'b');
        end


        clc
        fprintf('size... %d/%d \n', j, num_density);
        fprintf('\t (img) progress... %d/%d \n', i, num_image);
    end
    %%
    addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions')
    figure; hold on
    shadePlot(spatail_freq, mean(log10(col_freq_power), 1), std(log10(col_freq_power), [], 1)/sqrt(num_image), 'k');
    shadePlot(spatail_freq_r, mean(log10(col_freq_power_r), 1), std(log10(col_freq_power_r), [], 1)/sqrt(num_image), 'r');
    shadePlot(spatail_freq_rc, mean(log10(col_freq_power_rc), 1), std(log10(col_freq_power_rc), [], 1)/sqrt(num_image), 'b');
    
    %%
    keyboard; 

    Data{j, 1} = cData;

    dData = nan(num_mov_image, num_mov_sample, num_sim_cell);
    % moving natural image test
    for i = 1:num_mov_image
        mid = moving_test_image_id(i);
        file_name = sprintf('ResizeImg%d_%d.mat', FileTable(mid, 1), FileTable(mid, 2));
        load(fullfile(img_file_folder, file_name), 'img');
        img = squeeze(mean(double(img), 3))/255;
        img = 2*(img-0.5);
        img = img';

        mov_img = zeros(size(img));
        img_w = size(mov_img, 1);
        img_h = size(mov_img, 2);
        hids = round((0.5*(img_h-win_h)+1):(0.5*(img_h+win_h)));
        wids = round((0.5*(img_w-win_w)+1):(0.5*(img_w+win_w)));
        for k = 1:num_mov_sample
            if moving_x(k) <0
                x_loc_inds = 1:(img_w+moving_x(k));
                x_slc_inds = (-moving_x(k)+1):img_w;
            elseif moving_x(k) >0
                x_loc_inds = (moving_x(k)+1):img_w;
                x_slc_inds = 1:(img_w-moving_x(k));
            else
                x_loc_inds = 1:img_w;
                x_slc_inds = 1:img_w;
            end
            assert(length(x_loc_inds) == length(x_slc_inds));

            if moving_y(k) <0
                y_loc_inds = 1:(img_h+moving_y(k));
                y_slc_inds = (-moving_y(k)+1):img_h;
            elseif moving_y(k) >0
                y_loc_inds = (moving_y(k)+1):img_h;
                y_slc_inds = 1:(img_h-moving_y(k));
            else
                y_loc_inds = 1:img_h;
                y_slc_inds = 1:img_h;
            end
            assert(length(y_loc_inds) == length(y_slc_inds));
            mov_img(x_loc_inds, y_loc_inds) = img(x_slc_inds, y_slc_inds);
            % figure(1);
            % imagesc(mov_img);
            % title(sprintf('x:%d, y:%d', moving_x(k), moving_y(k)));
            % keyboard;

            mov_img = mov_img(wids, hids);

            cimg = mov_img(:)'*rf_maps;
            dData(i, k, :) = cimg;
            clc
            fprintf('size... %d/%d \n', j, num_density);
            fprintf('\t (mov) progress... %d/%d', i, num_mov_image);
        end
    end
    Data{j, 2} = dData;
end
%%
save_file_name = sprintf('%s_image_%s.mat', exp_cell_name, sim_data_time);
save_file_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\Simulation\Nat_Imgs';
save(fullfile(save_file_folder, save_file_name));
disp('Saved.')
%%
close all
plot_data = nan(1, num_mov_sample);
plot_data_s = nan(1, num_mov_sample);
plot_data_var = nan(num_density, num_mov_image);
explained_data = nan(num_density, 12);
explained_data_s = nan(num_density, 12);
n_noise_add = 5;
noise_level = 0.5;
sfscore = nan(size(score, 1), 2, n_noise_add, num_density);
avg_dist_data = nan(num_density, 1);
for j = 1:num_density
    size_id = j;
    cData = Data{size_id, 1};
    dData = Data{size_id, 2};
    points = Data{size_id, 3};
    mcData = mean(cData, 1);
    scData = std(cData, [], 1);
    % cData = cData-mcData;
    cData = (cData-mcData)./scData;
    [coeff,score,latent,tsquared,explained,mu] = pca(cData);
    for i = 1:n_noise_add
        noise_Data = cData+noise_level*randn(size(cData));
        tscore = noise_Data*coeff(:, 1:2);
        %noise_Data = (noise_Data-mean(noise_Data, 1))./std(noise_Data, [], 1);
        %[~,tscore] = pca(noise_Data);
        sfscore(:, :, i, j) = tscore(:, 1:2);
    end
    [explained_s, coeff_s, score_s] = spatialPCA(cData, points);
    explained_data(j, :) = explained(1:12);
    explained_data_s(j, :) = explained_s(1:12);

    dist_score = nan(num_mov_image, num_mov_sample);
    dist_score_s = nan(num_mov_image, num_mov_sample);
    % tdData = dData - mean(dData, [1 2]);
    % tdData = (dData-reshape(mcData, 1, 1, []));
    tdData = (dData-reshape(mcData, 1, 1, []))./reshape(scData, 1, 1, []);
    for i = 1:num_mov_image
        sdData = squeeze(tdData(i, :, :));
        d_score = sdData*coeff(:, 1:2);
        d_score_s = sdData*coeff_s(:, 1:2);
        dd_score = d_score - score(moving_test_image_id(i), 1:2);
        dd_score_s = d_score_s - score_s(moving_test_image_id(i), 1:2);
        dist_score(i, :) = sqrt(sum(dd_score.^2, 2));
        dist_score_s(i, :) = sqrt(sum(dd_score_s.^2, 2));
    end
    dist_mov = sqrt(moving_x.^2 + moving_y.^2);
    [~, sids] = sort(dist_mov);
    dist_score = dist_score(:, sids);
    dist_score_s = dist_score_s(:, sids);
    m_dist_score = mean(dist_score, 1);
    m_dist_score_s = mean(dist_score_s, 1);
    plot_data_var(j, :) = std(diff(dist_score, [], 2), [], 2);
    %
    dist_m = sqrt((score(:, 1)'-score(:, 1)).^2 + (score(:, 2)'-score(:, 2)).^2);
    gids = triu(ones(size(dist_m)), 1);
    avg_dist_m = mean(dist_m(gids == 1), 'all');
    avg_dist_data(j) = avg_dist_m;
    plot_data(j, :) = m_dist_score;
    plot_data_s(j, :) = m_dist_score_s;
    % plot_data(j, :) = m_dist_score/avg_dist_m;
    %
    figure;
    subplot(2, 2, 1)
    scatter(score(:, 1), score(:, 2), 5, 'k', 'filled');
    xlabel('PC1');
    ylabel('PC2');
    subplot(2, 2, 2); hold on
    plot(dist_mov(sids)*pix2um, dist_score');
    plot(dist_mov(sids)*pix2um, m_dist_score, 'k', 'LineWidth',2);
    xlabel('Moving distance (um)');
    ylabel('PC distance (12)');
    subplot(2, 2, 3)
    scatter(score_s(:, 1), score_s(:, 2), 5, 'k', 'filled');
    xlabel('PC1');
    ylabel('PC2');
    subplot(2, 2, 4); hold on
    plot(dist_mov(sids)*pix2um, dist_score_s');
    plot(dist_mov(sids)*pix2um, m_dist_score_s, 'k', 'LineWidth',2);
    xlabel('Moving distance (um)');
    ylabel('PC distance (12)');
end
%%
colors = parula(num_density);
clear density_legend
figure; hold on
for i = 1:num_density
    plot(1:12, explained_data(i, :), 'Color', colors(i, :));
    density_legend{i} = sprintf('%1G', target_density_list(i, 2));
end
xlabel('PC component ids');
ylabel('variance explained (%)');
legend(density_legend)
%%
pids = nchoosek(1:n_noise_add, 2);
num_pid = size(pids, 1);
dist_noise_add = nan(num_density, num_pid);
sfscore = real(sfscore);
for i = 1:num_density
    for j = 1:num_pid
        dist_noise_add(i, j) = squeeze(mean(sqrt(sum(sfscore(:, :, pids(j, 1), i)-sfscore(:, :, pids(j, 2), i).^2, 2)), 1));
    end

end
dist_noise_add = real(dist_noise_add)./avg_dist_data;
figure; hold on
clear tick_label
for i = 1:num_density
    scatter(i+0.7*(rand(1, num_pid)-0.5), dist_noise_add(i, :), 5, 'k', 'filled');
    plot(i+0.7*[-1 1]*0.5, median(dist_noise_add(i, :))*ones(1, 2), 'k');
    tick_label{i} = sprintf('%1G', target_density_list(i, 2));
end
xticks(1:num_density);
xticklabels(tick_label);
xlabel('RGC density (cell/mm2)');
ylabel('Noise influence (norm.)');
title(sprintf('Noise level: %.2G', noise_level))
    
%%
colors = parula(num_density);
close all
figure; 

for i = 1:num_density
    subplot(2, 2, 1);hold on
    plot(dist_mov(sids), plot_data(i, :), 'Color', colors(i, :));
    xlabel('Moved distance (um)');
    ylabel('Low dimension distance (PC unit)');
    legend(density_legend)
    subplot(2, 2, 2); hold on
    plot(dist_mov(sids), plot_data(i, :)./avg_dist_data(i), 'Color', colors(i, :));
    xlabel('Moved distance (um)');
    ylabel('Low dimension distance (norm.)');
    % scatter(i+0.7*(rand(1, num_mov_image)-0.5), plot_data_var(i, :), 5, 'k', 'filled');
    % plot(i+0.7*[-1 1]*0.5, median(plot_data_var(i, :))*ones(1, 2), 'k')
    subplot(2, 2, 3);hold on
    plot(dist_mov(sids), plot_data_s(i, :), 'Color', colors(i, :));
    xlabel('Moved distance (um)');
    ylabel('Low dimension distance (PC unit)');
    legend(density_legend)
    subplot(2, 2, 4);hold on
    plot(dist_mov(sids), plot_data_s(i, :)./avg_dist_data(i), 'Color', colors(i, :));
    xlabel('Moved distance (um)');
    ylabel('Low dimension distance (norm.)');
    % subplot(2, 2, 4); hold on
    % scatter(i+0.7*(rand(1, num_mov_image)-0.5), plot_data_var(i, :), 5, 'k', 'filled');
    % plot(i+0.7*[-1 1]*0.5, median(plot_data_var(i, :))*ones(1, 2), 'k')

end
%%
pids = nchoosek(1:5, 2);
keyboard;
%% preprocessing
% get nonlinear function (contrast -> firing rate)
[nl_fuc, divider]  = getNonlinearFunc(PBs, FRs);
% (1) find the center of RF (2) resample and (3) cut to a smaller cube
% extended size is determined by the size of position and RF
