% Parameters
show_cell_type = [0];        % Set to 0 or 1 to filter, or [] to ignore
show_location_type = [0];    % Set to 0 or 1 to filter, or [] to ignore
compare_type = '';           % Options: '', 'cell', 'location'
x_sample_range = -4:0.1:4;   % Should match the range used in NL_curves

numEntries = size(NL_curves, 1);
curves = {};
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
    curves{end+1} = NL_curves(i, :);
    titles{end+1} = sprintf('%s | %s | %s', data_sets{i}, cell_type{i}, location{i});
    filtered_indices(end+1) = i;
end

% Plotting
if strcmp(compare_type, 'cell') || strcmp(compare_type, 'location')
    leftCurves = {};
    rightCurves = {};
    leftTitles = {};
    rightTitles = {};
    for idx = 1:length(curves)
        orig_i = filtered_indices(idx);
        if strcmp(compare_type, 'cell')
            if cell_type_numeric(orig_i) == 0
                leftCurves{end+1} = curves{idx};
                leftTitles{end+1} = titles{idx};
            else
                rightCurves{end+1} = curves{idx};
                rightTitles{end+1} = titles{idx};
            end
        else
            if location_type_numeric(orig_i) == 0
                leftCurves{end+1} = curves{idx};
                leftTitles{end+1} = titles{idx};
            else
                rightCurves{end+1} = curves{idx};
                rightTitles{end+1} = titles{idx};
            end
        end
    end
    figure;
    subplot(1,2,1);
    hold on;
    for i = 1:length(leftCurves)
        plot(x_sample_range, leftCurves{i}, 'Color', [0.7 0.7 1]);
    end
    if ~isempty(leftCurves)
        avg_left = mean(cell2mat(leftCurves'),1);
        plot(x_sample_range, avg_left, 'b', 'LineWidth', 2);
    end
    hold off;
    title(sprintf('Type 0 (%s)', compare_type));
    xlabel('Generator signal'); ylabel('Firing rate');
    subplot(1,2,2);
    hold on;
    for i = 1:length(rightCurves)
        plot(x_sample_range, rightCurves{i}, 'Color', [1 0.7 0.7]);
    end
    if ~isempty(rightCurves)
        avg_right = mean(cell2mat(rightCurves'),1);
        plot(x_sample_range, avg_right, 'r', 'LineWidth', 2);
    end
    hold off;
    title(sprintf('Type 1 (%s)', compare_type));
    xlabel('Generator signal'); ylabel('Firing rate');
    sgtitle(sprintf('Nonlinearity Comparison by %s Type', compare_type));
else
    figure;
    hold on;
    for i = 1:length(curves)
        plot(x_sample_range, curves{i}, 'Color', [0.7 0.7 0.7]);
    end
    if ~isempty(curves)
        avg_all = mean(cell2mat(curves'),1);
        plot(x_sample_range, avg_all, 'k', 'LineWidth', 2);
    end
    hold off;
    title('Nonlinearity Curves');
    xlabel('Generator signal'); ylabel('Firing rate');
    legend([titles, {'Average'}], 'Interpreter', 'none');
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
save_file_name = fullfile(save_folder, sprintf('Nonlinearity%s%s_%s', filter_str, comp_str, datestr(now,'yyyymmdd_HHMMSS')));
print(gcf, save_file_name, '-depsc', '-painters');
print(gcf, save_file_name, '-dpng', '-r300');