% Parameters
zoom_in = 4;
show_cell_type = [1];        % Set to 0 or 1 to filter, or [] to ignore
show_location_type = [];    % Set to 0 or 1 to filter, or [] to ignore
compare_type = 'location';          % Options: '', 'cell', 'location'
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
    maxLen = max(length(leftImages), length(rightImages));
    figure;
    for i = 1:maxLen
        % Left column
        if i <= length(leftImages)
            subplot(maxLen, 2, (i - 1) * 2 + 1);
            if use_same_scale
                imagesc(leftImages{i}, [minVal, maxVal]);
            else
                imagesc(leftImages{i});
            end
            axis image off;
            title(leftTitles{i});
        end
        % Right column
        if i <= length(rightImages)
            subplot(maxLen, 2, (i - 1) * 2 + 2);
            if use_same_scale
                imagesc(rightImages{i}, [minVal, maxVal]);
            else
                imagesc(rightImages{i});
            end
            axis image off;
            title(rightTitles{i});
        end
    end
    
    sgtitle(sprintf('Zoomed x%d Comparison by %s Type', zoom_in, compare_type));
    colormap('jet');
    colorbar;

else
    % No comparison — regular grid layout
    numImages = length(images);
    cols = ceil(sqrt(numImages));
    rows = ceil(numImages / cols);
    
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
    end
    
    sgtitle(sprintf('Zoomed x%d Image Visualization', zoom_in));
    colormap('jet');
    colorbar;
end
