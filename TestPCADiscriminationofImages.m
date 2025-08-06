clear close all
%%
array_size = [900 1200];
num_of_points = 100;
x_bound = [100 500];
y_bound = [200 300];
random_seed = 42;
[x_indices, y_indices] = sample_indices_2d(array_size, num_of_points,...
     'XBound', x_bound, 'YBound', y_bound, 'RandSeed', random_seed);
figure; 
scatter(x_indices, y_indices, 5, 'b', 'filled');
%%
test_size = [10 20 40 80 120]; % 80 120];
extended_length = 100;
random_seed = 42;
x_bound = [extended_length+1 array_size(1)-extended_length];
y_bound = [extended_length+1 array_size(2)-extended_length];

%%
close all
% Specify the Excel file path
file_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\MEA\McGillDataset\MatImages';
file_name = 'FileTable.mat';
load(fullfile(file_folder, file_name), 'FileTable');
num_image = size(FileTable, 1);
rng(random_seed)
num_mov_image = 20;
moving_test_image_id = randperm(num_image, num_mov_image);
rng(random_seed)
num_mov_sample = 50;
moving_x = randsample(-99:99, num_mov_sample);
moving_x = sort(moving_x);
moving_y = randsample(-99:99, num_mov_sample);
moving_y = sort(moving_y);
clear Data
for j = 1:length(test_size)
    cData = nan(num_image, test_size(j));
    [x_indices, y_indices] = sample_indices_2d(array_size, test_size(j),...
        'XBound', x_bound, 'YBound', y_bound, 'RandSeed', random_seed);
    indices = sub2ind(array_size, x_indices, y_indices);
    % every natural image
    for i = 1:num_image
        file_name = sprintf('ResizeImg%d_%d.mat', FileTable(i, 1), FileTable(i, 2));
        load(fullfile(file_folder, file_name), 'img');
        img = squeeze(mean(double(img), 3))/255;
        % figure(1);
        % imagesc(img, [0 1]);
        cData(i, :) = img(indices);
        clc
        fprintf('size... %d/%d \n', j, length(test_size));
        fprintf('\t (img) progress... %d/%d \n', i, num_image);
    end
    
    Data{j, 1} = cData;
    dData = nan(num_mov_image, num_mov_sample, test_size(j));
    % moving natural image test
    for i = 1:num_mov_image
        mid = moving_test_image_id(i);
        file_name = sprintf('ResizeImg%d_%d.mat', FileTable(mid, 1), FileTable(mid, 2));
        load(fullfile(file_folder, file_name), 'img');
        img = squeeze(mean(double(img), 3))/255;
        mov_img = zeros(size(img));
        for k = 1:num_mov_sample
            if moving_x(k) <0
                x_loc_inds = 1:(array_size(1)+moving_x(k));
                x_slc_inds = (-moving_x(k)+1):array_size(1);
            elseif moving_x(k) >0
                x_loc_inds = (moving_x(k)+1):array_size(1);
                x_slc_inds = 1:(array_size(1)-moving_x(k));
            else
                x_loc_inds = 1:array_size(1);
                x_slc_inds = 1:array_size(1);
            end
            assert(length(x_loc_inds) == length(x_slc_inds));

            if moving_y(k) <0
                y_loc_inds = 1:(array_size(2)+moving_y(k));
                y_slc_inds = (-moving_y(k)+1):array_size(2);
            elseif moving_y(k) >0
                y_loc_inds = (moving_y(k)+1):array_size(2);
                y_slc_inds = 1:(array_size(2)-moving_y(k));
            else
                y_loc_inds = 1:array_size(2);
                y_slc_inds = 1:array_size(2);
            end
            assert(length(y_loc_inds) == length(y_slc_inds));
            mov_img(x_loc_inds, y_loc_inds) = img(x_slc_inds, y_slc_inds);
            % figure(1);
            % imagesc(mov_img);
            % title(sprintf('x:%d, y:%d', moving_x(k), moving_y(k)));
            % keyboard;

            dData(i, k, :) = mov_img(indices);
            clc
            fprintf('size... %d/%d \n', j, length(test_size));
            fprintf('\t (mov) progress... %d/%d', i, num_mov_image);
        end
    end
    Data{j, 2} = dData;
end

%%
save_file_name = '0820202402_image_pixel_test.mat';
save_file_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\Simulation\Toy_nat_imgs';
save(fullfile(save_file_folder, save_file_name));
disp('Saved.')
%%
close all
plot_data = nan(length(test_size), num_mov_sample);
for j = 1:length(test_size)
    size_id = j;
    cData = Data{size_id, 1};
    dData = Data{size_id, 2};
    [coeff,score,latent,tsquared,explained,mu] = pca(cData);
    
    dist_score = nan(num_mov_image, num_mov_sample);
    tdData = dData - mean(dData, [1 2]);
    for i = 1:num_mov_image
        sdData = squeeze(tdData(i, :, :));
        d_score = sdData*coeff(:, 1:2);
        dd_score = d_score - score(moving_test_image_id(i), 1:2);
        dist_score(i, :) = sqrt(sum(dd_score.^2, 2));
    end
    dist_mov = sqrt(moving_x.^2 + moving_y.^2);
    [~, sids] = sort(dist_mov);
    dist_score = dist_score(:, sids);
    m_dist_score = mean(dist_score, 1);

    %
    dist_m = sqrt((score(:, 1)'-score(:, 1)).^2 + (score(:, 2)'-score(:, 2)).^2);
    gids = triu(ones(size(dist_m)), 1);
    avg_dist_m = mean(dist_m(gids == 1), 'all');
    plot_data(j, :) = m_dist_score/avg_dist_m;
    %
    figure;
    subplot(1, 2, 1)
    scatter(score(:, 1), score(:, 2), 5, 'k', 'filled');
    subplot(1, 2, 2); hold on
    plot(dist_mov(sids), dist_score');
    plot(dist_mov(sids), m_dist_score, 'k', 'LineWidth',2);
end

colors = parula(length(test_size));
figure; hold on
for i = 1:length(test_size)
    plot(dist_mov(sids), plot_data(i, :), 'Color', colors(i, :));
end