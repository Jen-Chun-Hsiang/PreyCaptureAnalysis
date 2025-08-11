% Parameters
show_cell_type = [1];        % Set to 0 or 1 to filter, or [] to ignore
show_location_type = [0];    % Set to 0 or 1 to filter, or [] to ignore
compare_type = '';           % Options: '', 'cell', 'location'
normalize_peak = true;      % Set to true to normalize each filter to its peak (abs)
sampling_freq = 100;         % Hz

numEntries = length(Data);
filters = {};
titles = {};
filtered_indices = [];

% Conflict resolution
if strcmp(compare_type, 'cell') && ~isempty(show_cell_type)
    warning('Ignoring show_cell_type due to compare_type = "cell".');
    show_cell_type = [];
end
if strcmp(compare_type, 'location') && ~isempty(show_location_type)
    warning('Ignoring show_location_type due to compare_type = "location".');
    show_location_type = [];
end

% Filtering and extraction
for i = 1:numEntries
    ct = cell_type_numeric(i);
    lt = location_type_numeric(i);
    if ~isempty(show_cell_type) && ct ~= show_cell_type
        continue;
    end
    if ~isempty(show_location_type) && lt ~= show_location_type
        continue;
    end
    tfilt = Data{i}.tRF(:)';
    if normalize_peak
        tfilt = tfilt / max(abs(tfilt));
    end
    filters{end+1} = tfilt;
    titles{end+1} = sprintf('%s | %s | %s', data_sets{i}, cell_type{i}, location{i});
    filtered_indices(end+1) = i;
end

% Time axis (backward in time)
if ~isempty(filters)
    nT = length(filters{1});
    t = -(nT-1)/sampling_freq : 1/sampling_freq : 0;
else
    t = [];
end

% Plotting
if strcmp(compare_type, 'cell') || strcmp(compare_type, 'location')
    leftFilters = {};
    rightFilters = {};
    leftTitles = {};
    rightTitles = {};
    for idx = 1:length(filters)
        orig_i = filtered_indices(idx);
        if strcmp(compare_type, 'cell')
            if cell_type_numeric(orig_i) == 0
                leftFilters{end+1} = filters{idx};
                leftTitles{end+1} = titles{idx};
            else
                rightFilters{end+1} = filters{idx};
                rightTitles{end+1} = titles{idx};
            end
        else
            if location_type_numeric(orig_i) == 0
                leftFilters{end+1} = filters{idx};
                leftTitles{end+1} = titles{idx};
            else
                rightFilters{end+1} = filters{idx};
                rightTitles{end+1} = titles{idx};
            end
        end
    end
    figure;
    subplot(1,2,1);
    hold on;
    for i = 1:length(leftFilters)
        plot(t, leftFilters{i}, 'Color', [0.7 0.7 1]);
    end
    if ~isempty(leftFilters)
        avg_left = mean(cell2mat(leftFilters'),1);
        plot(t, avg_left, 'b', 'LineWidth', 2);
    end
    hold off;
    title(sprintf('Type 0 (%s)', compare_type));
    xlabel('Time before spike (s)'); ylabel('Filter');
    if normalize_peak
        ylim([-1 1]);
        yticks(-1:0.5:1);
        yticklabels({'-1','', '0', '', '1'});
        xlim([-0.5 0]);
        xticks([-0.5 -0.25 0]);
        xticklabels({'-0.5','-0.25','0'});
    end
    subplot(1,2,2);
    hold on;
    for i = 1:length(rightFilters)
        plot(t, rightFilters{i}, 'Color', [1 0.7 0.7]);
    end
    if ~isempty(rightFilters)
        avg_right = mean(cell2mat(rightFilters'),1);
        plot(t, avg_right, 'r', 'LineWidth', 2);
    end
    hold off;
    title(sprintf('Type 1 (%s)', compare_type));
    xlabel('Time before spike (s)'); ylabel('Filter');
    if normalize_peak
        ylim([-1 1]);
        yticks(-1:0.5:1);
        yticklabels({'-1','', '0', '', '1'});
        xlim([-0.5 0]);
        xticks([-0.5 -0.25 0]);
        xticklabels({'-0.5','-0.25','0'});
    end
    sgtitle(sprintf('Temporal Filter Comparison by %s Type', compare_type));
else
    figure;
    hold on;
    for i = 1:length(filters)
        plot(t, filters{i}, 'Color', [0.7 0.7 0.7]);
    end
    if ~isempty(filters)
        avg_all = mean(cell2mat(filters'),1);
        plot(t, avg_all, 'k', 'LineWidth', 2);
    end
    hold off;
    title('Temporal Filters');
    xlabel('Time before spike (s)'); ylabel('Filter');
    legend([titles, {'Average'}], 'Interpreter', 'none');
    if normalize_peak
        ylim([-1 1]);
        yticks(-1:0.5:1);
        yticklabels({'-1','', '0', '', '1'});
        xlim([-0.5 0]);
        xticks([-0.5 -0.25 0]);
        xticklabels({'-0.5','-0.25','0'});
    end
end

% Save figure
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
norm_str = '';
if normalize_peak
    norm_str = '_normpeak';
end

save_file_name = fullfile(save_folder, sprintf('TemporalFilters%s%s%s_%s', filter_str, comp_str, norm_str, datestr(now,'yyyymmdd_HHMMSS')));
print(gcf, save_file_name, '-depsc', '-painters');
print(gcf, save_file_name, '-dpng', '-r300');
