% Script: FindThreshold_MeanMinusKStd.m
% Purpose: For each image, use its fitted Gaussian parameters (gauss_est) and image (Data{k}.stdSTA')
% to compute a mask, then find a robust threshold using mean minus k*std within the mask.
% The threshold is then reported for each image.

% --- User settings ---
std_plus = 2; % Number of standard deviations below the mean
num_gauss = 1; % Number of Gaussians used in fitting

% --- Load processed results ---

% processedFile = 'GaussianFitting_processed_080725_1.mat';
% load(processedFile, 'gauss_est');

% folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
num_set = size(gauss_est, 1);
thresholds = nan(num_set, 1);
rf_pixels = nan(num_set, 1); % To store RF sizes

for k = 1:num_set
    % file_name = sprintf('%s.mat', data_sets{k});
    % S = load(fullfile(folder_name, file_name), 'stdSTA');
    image = Data{k}.stdSTA';
    params = gauss_est(k, :);
    % Create mask from Gaussian fit (ellipse at 2*sigma)
    cx = params(1); cy = params(2); sx = params(3); sy = params(4); theta = params(5);
    [X, Y] = meshgrid(1:size(image, 2), 1:size(image, 1));
    Xc = X - cx; Yc = Y - cy;
    Xr = Xc * cos(theta) + Yc * sin(theta);
    Yr = -Xc * sin(theta) + Yc * cos(theta);
    mask = (Xr.^2/(2*sx^2) + Yr.^2/(2*sy^2)) <= 2^2; % 2*sigma ellipse
    % Find threshold within mask
    mask_pixels = image(~mask);
    mu = mean(mask_pixels);
    sigma = std(mask_pixels);
    threshold = mu + std_plus*sigma;
    thresholds(k) = threshold;
    fprintf('Image %d: threshold = %.4f\n', k, threshold);
    % Binary mask
    bw = image > threshold;
    % Find connected components
    CC = bwconncomp(bw);
    stats = regionprops(CC, 'Area');
    if isempty(stats)
        largest_area = 0;
    else
        areas = [stats.Area];
        [largest_area, idx] = max(areas);
        % Optionally, get the mask for the largest group
        largest_mask = false(size(bw));
        largest_mask(CC.PixelIdxList{idx}) = true;
    end
    rf_pixels(k) = largest_area;
    fprintf('Image %d: largest connected group size = %d pixels\n', k, largest_area);
    % figure; imagesc(largest_mask); title(sprintf('Largest group, %d pixels', largest_area));
    % keyboard
end
rf_pixels*4.375^2
mean(rf_pixels*4.375^2)
mean(sqrt(rf_pixels*4.375^2*4/pi))
processedFile = fullfile(folder_name, process_version);
if exist(processedFile, 'file')
    save(processedFile, 'rf_pixels', '-append');
end

% Save thresholds
% save('Thresholds_MeanMinusKStd.mat', 'thresholds', 'k');
