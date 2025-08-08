% Parameters
zoom_in = 4;
show_cell_type = [0];        % Set to 0 or 1 to filter, or [] to ignore
show_location_type = [0];    % Set to 0 or 1 to filter, or [] to ignore
compare_type = '';          % Options: '', 'cell', 'location'
use_same_scale = false;      % Set to false to scale each image individually

numEntries = length(Data);

images = {};
titles = {};
filtered_indices = [];
minVal = inf;
maxVal = -inf;

% Conflict resolution
if strcmp(compare_type, 'cell') && ~isempty(show_cell_type)
    warning('Ignoring show_cell_type due to compare_type = "cell".');
    show_cell_type = [];
end
if strcmp(compare_type, 'location') && ~isempty(show_location_type)
    warning('Ignoring show_location_type due to compare_type = "location".');
    show_location_type = [];
end

% Filtering and image extraction
for i = 1:numEntries
    ct = cell_type_numeric(i);
    lt = location_type_numeric(i);
    
    % Apply filters
    if ~isempty(show_cell_type) && ct ~= show_cell_type
        continue;
    end
    if ~isempty(show_location_type) && lt ~= show_location_type
        continue;
    end
    
    img = Data{i}.stdSTA';
    [h, w] = size(img);
    zh = round(h / zoom_in / 2);
    zw = round(w / zoom_in / 2);
    center_h = round(h / 2);
    center_w = round(w / 2);
    cropped = img(center_h - zh + 1 : center_h + zh, ...
                  center_w - zw + 1 : center_w + zw);
    
    images{end+1} = cropped;
    titles{end+1} = sprintf('%s | %s | %s', data_sets{i}, cell_type{i}, location{i});
    filtered_indices(end+1) = i;
    minVal = min(minVal, min(cropped(:)));
    maxVal = max(maxVal, max(cropped(:)));
end

% Comparison grouping
if strcmp(compare_type, 'cell') || strcmp(compare_type, 'location')
    leftImages = {};
    rightImages = {};
    leftTitles = {};
    rightTitles = {};
    for idx = 1:length(images)
        orig_i = filtered_indices(idx);
        if strcmp(compare_type, 'cell')
            if cell_type_numeric(orig_i) == 0
                leftImages{end+1} = images{idx};
                leftTitles{end+1} = titles{idx};
            else
                rightImages{end+1} = images{idx};
                rightTitles{end+1} = titles{idx};
            end
        else
            if location_type_numeric(orig_i) == 0
                leftImages{end+1} = images{idx};
                leftTitles{end+1} = titles{idx};
            else
                rightImages{end+1} = images{idx};
                rightTitles{end+1} = titles{idx};
            end
        end
    end
    
    % Determine layout
    leftCols = 2;
    rightCols = 2;
    totalCols = leftCols + rightCols; % 4 columns
    maxRows = max(ceil(length(leftImages)/leftCols), ceil(length(rightImages)/rightCols));
    figure;
    % Plot left images (columns 1-2)
    for i = 1:length(leftImages)
        row = ceil(i/leftCols);
        col = mod(i-1, leftCols) + 1;
        subplot(maxRows, totalCols, (row-1)*totalCols + col);
        if use_same_scale
            imagesc(leftImages{i}, [minVal, maxVal]);
        else
            imagesc(leftImages{i});
        end
        axis image off;
        title(leftTitles{i});
        % Set colormap based on cell_type_numeric
        orig_i = filtered_indices(i);
        if cell_type_numeric(orig_i) == 1
            colormap(gray);
        else
            colormap(flipud(gray));
        end
        % Overlay ellipse from gauss_est (center aligned to cropped image)
        if exist('gauss_est', 'var') && size(gauss_est, 1) >= orig_i
            % Get full image size and cropping info
            img_full = Data{orig_i}.stdSTA';
            [h, w] = size(img_full);
            zh = round(h / zoom_in / 2);
            zw = round(w / zoom_in / 2);
            center_h = round(h / 2);
            center_w = round(w / 2);
            crop_center_h = zh + 1;
            crop_center_w = zw + 1;
            % Offset between full image center and cropped image center
            offset_x = crop_center_w - center_w;
            offset_y = crop_center_h - center_h;
            % Ellipse parameters
            x_mean = gauss_est(orig_i, 1);
            y_mean = gauss_est(orig_i, 2);
            sigma_x = gauss_est(orig_i, 3);
            sigma_y = gauss_est(orig_i, 4);
            theta = gauss_est(orig_i, 5);
            t_ellipse = linspace(0, 2*pi, 100);
            x_ellipse = sigma_x * cos(t_ellipse);
            y_ellipse = sigma_y * sin(t_ellipse);
            R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
            ellipse_coords = R * [x_ellipse; y_ellipse];
            % Shift ellipse center to cropped image center
            x_plot = ellipse_coords(1, :) + (x_mean - center_w + crop_center_w);
            y_plot = ellipse_coords(2, :) + (y_mean - center_h + crop_center_h);
            hold on;
            plot(x_plot, y_plot, 'b', 'LineWidth', 2);
            hold off;
        end
    end
    % Plot right images (columns 3-4)
    for i = 1:length(rightImages)
        row = ceil(i/rightCols);
        col = mod(i-1, rightCols) + 1 + leftCols; % shift to columns 3-4
        subplot(maxRows, totalCols, (row-1)*totalCols + col);
        if use_same_scale
            imagesc(rightImages{i}, [minVal, maxVal]);
        else
            imagesc(rightImages{i});
        end
        axis image off;
        title(rightTitles{i});
        % Set colormap based on cell_type_numeric
        orig_i = filtered_indices(i + length(leftImages));
        if cell_type_numeric(orig_i) == 1
            colormap(gray);
        else
            colormap(flipud(gray));
        end
        % Overlay ellipse from gauss_est (center aligned to cropped image)
        if exist('gauss_est', 'var') && size(gauss_est, 1) >= orig_i
            % Get full image size and cropping info
            img_full = Data{orig_i}.stdSTA';
            [h, w] = size(img_full);
            zh = round(h / zoom_in / 2);
            zw = round(w / zoom_in / 2);
            center_h = round(h / 2);
            center_w = round(w / 2);
            crop_center_h = zh + 1;
            crop_center_w = zw + 1;
            % Offset between full image center and cropped image center
            offset_x = crop_center_w - center_w;
            offset_y = crop_center_h - center_h;
            % Ellipse parameters
            x_mean = gauss_est(orig_i, 1);
            y_mean = gauss_est(orig_i, 2);
            sigma_x = gauss_est(orig_i, 3);
            sigma_y = gauss_est(orig_i, 4);
            theta = gauss_est(orig_i, 5);
            t_ellipse = linspace(0, 2*pi, 100);
            x_ellipse = sigma_x * cos(t_ellipse);
            y_ellipse = sigma_y * sin(t_ellipse);
            R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
            ellipse_coords = R * [x_ellipse; y_ellipse];
            % Shift ellipse center to cropped image center
            x_plot = ellipse_coords(1, :) + (x_mean - center_w + crop_center_w);
            y_plot = ellipse_coords(2, :) + (y_mean - center_h + crop_center_h);
            hold on;
            plot(x_plot, y_plot, 'b', 'LineWidth', 2);
            hold off;
        end
    end
    sgtitle(sprintf('Zoomed x%d Comparison by %s Type', zoom_in, compare_type));
    % colorbar;

else
    % No comparison — regular grid layout
    numImages = length(images);
    cols = ceil(sqrt(numImages+1)); % +1 for the extra image
    rows = ceil((numImages+1) / cols);
    if numImages == rows*cols
        rows = rows + 1;
    end
    figure;
    for i = 1:numImages
        subplot(rows, cols, i);
        if use_same_scale
            imagesc(images{i}, [minVal, maxVal]);
        else
            imagesc(images{i});
        end
        axis image off;
        title(titles{i});
        % Set colormap based on cell_type_numeric
        orig_i = filtered_indices(i);
        if cell_type_numeric(orig_i) == 1
            colormap(gray);
        else
            colormap(flipud(gray));
        end
        % Overlay ellipse from gauss_est (center aligned to cropped image)
        if exist('gauss_est', 'var') && size(gauss_est, 1) >= orig_i
            img_full = Data{orig_i}.stdSTA';
            [h, w] = size(img_full);
            zh = round(h / zoom_in / 2);
            zw = round(w / zoom_in / 2);
            center_h = round(h / 2);
            center_w = round(w / 2);
            crop_center_h = zh + 1;
            crop_center_w = zw + 1;
            x_mean = gauss_est(orig_i, 1);
            y_mean = gauss_est(orig_i, 2);
            sigma_x = gauss_est(orig_i, 3);
            sigma_y = gauss_est(orig_i, 4);
            theta = gauss_est(orig_i, 5);
            t_ellipse = linspace(0, 2*pi, 100);
            x_ellipse = 2*sigma_x * cos(t_ellipse);
            y_ellipse = 2*sigma_y * sin(t_ellipse);
            R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
            ellipse_coords = R * [x_ellipse; y_ellipse];
            x_plot = ellipse_coords(1, :) + (x_mean - center_w + crop_center_w);
            y_plot = ellipse_coords(2, :) + (y_mean - center_h + crop_center_h);
            hold on;
            plot(x_plot, y_plot, 'b', 'LineWidth', 2);
            hold off;
        end
    end
    % Add one more image (duplicate the first one) and plot all ellipses of specified cell type and location
    subplot(rows, cols, numImages+1);
    imagesc(images{1});
    axis image off;
    title('All ellipses overlay');
    orig_i_first = filtered_indices(1);
    if cell_type_numeric(orig_i_first) == 1
        colormap(gray);
    else
        colormap(flipud(gray));
    end
    hold on;
    for idx = 1:length(filtered_indices)
        orig_i = filtered_indices(idx);
        img_full = Data{orig_i}.stdSTA';
        [h, w] = size(img_full);
        zh = round(h / zoom_in / 2);
        zw = round(w / zoom_in / 2);
        center_h = round(h / 2);
        center_w = round(w / 2);
        crop_center_h = zh + 1;
        crop_center_w = zw + 1;
        x_mean = gauss_est(orig_i, 1);
        y_mean = gauss_est(orig_i, 2);
        sigma_x = gauss_est(orig_i, 3);
        sigma_y = gauss_est(orig_i, 4);
        theta = gauss_est(orig_i, 5);
        t_ellipse = linspace(0, 2*pi, 100);
        x_ellipse = 2*sigma_x * cos(t_ellipse);
        y_ellipse = 2*sigma_y * sin(t_ellipse);
        R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        ellipse_coords = R * [x_ellipse; y_ellipse];
        x_plot = ellipse_coords(1, :) + (x_mean - center_w + crop_center_w);
        y_plot = ellipse_coords(2, :) + (y_mean - center_h + crop_center_h);
        plot(x_plot, y_plot, 'b', 'LineWidth', 2);
    end
    hold off;
    sgtitle(sprintf('Zoomed x%d Image Visualization', zoom_in));
    % colorbar;
end
%%
% Build descriptive filename for figure
filter_str = '';
if ~isempty(show_cell_type)
    filter_str = sprintf('_celltype%d', show_cell_type);
end
if ~isempty(show_location_type)
    filter_str = [filter_str, sprintf('_locationtype%d', show_location_type)];
end
comp_str = '';
if ~isempty(compare_type)
    comp_str = ['_', compare_type];
end
save_file_name = fullfile(save_folder, sprintf('RF_Images_Zoomx%d%s%s_%s', zoom_in, filter_str, comp_str, datestr(now,'yyyymmdd_HHMMSS')));
print(gcf, save_file_name, '-depsc', '-painters'); % EPS format
print(gcf, save_file_name, '-dpng', '-r300'); % PNG, 600 dpi


%%
% Script to plot ellipses for each Gaussian fit in gauss_est
% Assumes gauss_est is already in the workspace

% figure; hold on; axis equal;
% n = size(gauss_est, 1);

% t = linspace(0, 2*pi, 100); % parameter for ellipse

% for i = 1:n
%     x_mean = gauss_est(i, 1);
%     y_mean = gauss_est(i, 2);
%     sigma_x = gauss_est(i, 3);
%     sigma_y = gauss_est(i, 4);
%     theta = gauss_est(i, 5);

%     % Ellipse before rotation
%     x_ellipse = sigma_x * cos(t);
%     y_ellipse = sigma_y * sin(t);

%     % Rotation matrix
%     R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
%     ellipse_coords = R * [x_ellipse; y_ellipse];

%     % Shift to center
%     x_plot = ellipse_coords(1, :) + x_mean;
%     y_plot = ellipse_coords(2, :) + y_mean;

%     plot(x_plot, y_plot, 'LineWidth', 2);
% end

% xlabel('X');
% ylabel('Y');
% title('Ellipses from gauss\_est');
% hold off;