% Simplified spot size analysis: Calculate Size Index (SI) for different groups
% SI = (S-C)/C where:
%   S = mean response of spot sizes 800, 1200
%   C = mean response of top two spot sizes < 800
%
% Requires: struct 'a' in workspace with fields:
%   - x1 (T x 1), y_labels (1 x S), and group arrays (T x (S * Ncells))
% Example groups: 'DN_ONSus_RF_UV', 'DN_ONSus_RF_GRN'

clc; close all;
test_type = 'OFF';
percent_peak = 0.85;

% Save figure folder
save_fig_folder = './Figures/SpotSizeSimple/';
if ~exist(save_fig_folder, 'dir')
    mkdir(save_fig_folder);
end
a = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\Spots\SusAlpha_RawData.mat');

% ---------------- Configuration ----------------
if ~exist('a','var')
    error('Struct ''a'' not found in workspace.');
end
switch test_type
    case 'ON'
        groups = {'AcuteZoneDT_ONSus_RF_GRN','DN_ONSus_RF_GRN'};  % edit/extend as needed
        stim_idx = 110:300;  % during-stimulus period
        % stim_idx = 110:150;  % during-stimulus period
    case 'OFF'  
        groups = {'AcuteZoneDT_OFFSus_RF_GRN','DN_OFFSus_RF_GRN'};  % edit/extend as needed
        stim_idx = 310:500;  % during-stimulus period
        % stim_idx = 310:350;  % during-stimulus period
end
sizes = a.y_labels(:)';                  % diameter labels
nSizes = numel(sizes);
T = numel(a.x1);

% Ensure sizes are ascending
[sorted_sizes, sort_idx] = sort(sizes, 'ascend');

% Define size categories for SI calculation
large_sizes = [800, 1200];  % Spot sizes for S (surround suppression sizes)
small_size_threshold = 800; % Threshold for C (center sizes)

fprintf('Size categories for SI calculation:\n');
fprintf('  Large sizes (S): %s\n', mat2str(large_sizes));
fprintf('  Small size threshold (C): < %d\n', small_size_threshold);
fprintf('  Available sizes: %s\n', mat2str(sorted_sizes));

% ---------------- Analysis ----------------
results = struct();
all_SI_values = [];
all_group_labels = [];
group_colors = lines(numel(groups));

for g = 1:numel(groups)
    gname = groups{g};
    if ~isfield(a, gname)
        warning('Group "%s" not found in struct a. Skipping.', gname);
        continue;
    end

    M = a.(gname);            % T x (nSizes * nCells)
    [Tg, nCol] = size(M);

    % Basic checks
    if Tg ~= T
        error('Time length mismatch for group %s (got %d, expected %d).', gname, Tg, T);
    end
    if mod(nCol, nSizes) ~= 0
        error('Columns (%d) not divisible by number of sizes (%d) for group %s.', nCol, nSizes, gname);
    end

    nCells = nCol / nSizes;
    fprintf('\nProcessing group %s: %d cells\n', gname, nCells);

    % Mean across time during stimulus for each column (cell x size)
    stim_mean_by_col = mean(M(stim_idx, :), 1, 'omitnan');  % 1 x (nSizes*nCells)

    % Reshape into [nSizes x nCells]; columns order is per cell then sizes
    S_matrix = reshape(stim_mean_by_col, [nSizes, nCells]);

    % Reorder sizes to be ascending (and reorder S_matrix accordingly)
    S_matrix = S_matrix(sort_idx, :);

    % Calculate SI index for each cell
    SI_values = nan(nCells, 1);
    
    for c = 1:nCells
        cell_responses = S_matrix(:, c);  % Response of cell c to different sizes
        
        % Find indices for large sizes (S)
        large_size_indices = ismember(sorted_sizes, large_sizes);
        if sum(large_size_indices) == 0
            fprintf('  Cell %d: No large sizes found, skipping\n', c);
            continue;
        end
        
        % Find indices for small sizes (C)
        small_size_indices = sorted_sizes < small_size_threshold;
        if sum(small_size_indices) < 2
            fprintf('  Cell %d: Less than 2 small sizes found, skipping\n', c);
            continue;
        end
        
        % Calculate S: mean response to large sizes
        S_responses = cell_responses(large_size_indices);
        S_mean = mean(S_responses, 'omitnan');
        
        % Calculate C: peak response and closest larger adjacent response from small sizes
        C_responses = cell_responses(small_size_indices);
        small_sizes_for_this_cell = sorted_sizes(small_size_indices);
        
        % Remove NaN values for peak finding
        valid_indices = ~isnan(C_responses);
        if sum(valid_indices) < 2
            fprintf('  Cell %d: Less than 2 valid small size responses, skipping\n', c);
            continue;
        end
        
        valid_C_responses = C_responses(valid_indices);
        valid_small_sizes = small_sizes_for_this_cell(valid_indices);
        
        % Find peak response
        [peak_value, peak_idx] = max(valid_C_responses);
        peak_size = valid_small_sizes(peak_idx);
        
        % Find adjacent spot sizes to the peak
        peak_size_original_idx = find(sorted_sizes == peak_size);
        
        % Look for adjacent sizes (both smaller and larger than peak)
        adjacent_responses = [];
        adjacent_info = {};
        
        % Check size immediately before peak (if exists and < small_size_threshold)
        if peak_size_original_idx > 1
            prev_size = sorted_sizes(peak_size_original_idx - 1);
            if prev_size < small_size_threshold && ~isnan(cell_responses(peak_size_original_idx - 1))
                adjacent_responses(end+1) = cell_responses(peak_size_original_idx - 1);
                adjacent_info{end+1} = sprintf('size %d (smaller)', prev_size);
            end
        end
        
        % Check size immediately after peak (if exists and < small_size_threshold)
        if peak_size_original_idx < length(sorted_sizes)
            next_size = sorted_sizes(peak_size_original_idx + 1);
            if next_size < small_size_threshold && ~isnan(cell_responses(peak_size_original_idx + 1))
                adjacent_responses(end+1) = cell_responses(peak_size_original_idx + 1);
                adjacent_info{end+1} = sprintf('size %d (larger)', next_size);
            end
        end
        
        % Calculate C_mean based on available adjacent responses
        if length(adjacent_responses) >= 1
            % Choose the largest adjacent response
            [largest_adjacent, largest_idx] = max(adjacent_responses);
            C_mean = (peak_value + largest_adjacent) / 2;
            
            % Debug info for first few cells
            if c <= 3
                fprintf('  Cell %d: Peak=%.2f at size %d, Adjacent=%.2f (%s), C_mean=%.2f\n', ...
                    c, peak_value, peak_size, largest_adjacent, adjacent_info{largest_idx}, C_mean);
            end
        else
            fprintf('  Cell %d: No valid adjacent responses to peak, skipping\n', c);
            continue;
        end
        
        % Calculate SI = (S-C)/C
        if C_mean > 0
            SI_values(c) = (S_mean - C_mean) / C_mean;
        else
            fprintf('  Cell %d: C_mean <= 0, skipping\n', c);
            continue;
        end
        
        % Debug output for first few cells
        if c <= 3
            fprintf('  Cell %d: S=%.2f, C=%.2f, SI=%.3f\n', c, S_mean, C_mean, SI_values(c));
        end
    end
    
    % Remove NaN values
    valid_SI = ~isnan(SI_values);
    valid_SI_values = SI_values(valid_SI);
    
    fprintf('  Valid SI calculations: %d/%d cells\n', sum(valid_SI), nCells);
    fprintf('  SI range: %.3f to %.3f\n', min(valid_SI_values), max(valid_SI_values));
    fprintf('  Mean SI: %.3f ± %.3f (SEM)\n', mean(valid_SI_values), std(valid_SI_values)/sqrt(length(valid_SI_values)));
    
    % Store results
    results.(gname).SI_values = SI_values;
    results.(gname).valid_SI = valid_SI;
    results.(gname).mean_SI = mean(valid_SI_values);
    results.(gname).sem_SI = std(valid_SI_values) / sqrt(length(valid_SI_values));
    results.(gname).std_SI = std(valid_SI_values);
    results.(gname).n_cells = length(valid_SI_values);
    results.(gname).sizes = sorted_sizes;
    results.(gname).responses_matrix = S_matrix;
    
    % Collect data for combined plot
    all_SI_values = [all_SI_values; valid_SI_values];
    all_group_labels = [all_group_labels; repmat(g, length(valid_SI_values), 1)];
end

% ---------------- Plotting ----------------
if ~isempty(all_SI_values)
    figure('Color', 'w', 'Position', [100, 100, 800, 600]);
    
    % Calculate bar positions
    group_positions = 1:numel(groups);
    bar_width = 0.6;
    
    % Prepare data for bar plot
    mean_values = nan(numel(groups), 1);
    sem_values = nan(numel(groups), 1);
    
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(results, gname)
            mean_values(g) = results.(gname).mean_SI;
            sem_values(g) = results.(gname).sem_SI;
        end
    end
    
    % Create bar plot with error bars
    hold on;
    bars = bar(group_positions, mean_values, bar_width, 'FaceAlpha', 0.7);
    
    % Add error bars
    errorbar(group_positions, mean_values, sem_values, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'CapSize', 8);
    
    % Add individual data points
    jitter_width = 0.15;
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(results, gname)
            valid_SI_values = results.(gname).SI_values(results.(gname).valid_SI);
            
            % Add jitter to x-positions
            n_points = length(valid_SI_values);
            x_jitter = (rand(n_points, 1) - 0.5) * jitter_width * 2;
            x_positions = repmat(group_positions(g), n_points, 1) + x_jitter;
            
            % Plot individual points
            scatter(x_positions, valid_SI_values, 40, group_colors(g,:), 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8);
        end
    end
    
    % Customize plot
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    xlabel('Cell Groups');
    ylabel('Size Index (SI)');
    title('Size Index: (Large Spot Response - Center Response) / Center Response');
    grid on;
    
    % Add horizontal line at SI = 0
    yline(0, 'k--', 'LineWidth', 1, 'Alpha', 0.5);
    
    % Add text annotation explaining SI calculation
    annotation_text = sprintf(['SI = (S - C) / C\n' ...
        'S = mean response to spots %s\n' ...
        'C = mean of peak response and largest adjacent response\n' ...
        'for spots < %d'], ...
        mat2str(large_sizes), small_size_threshold);
    
    text(0.02, 0.98, annotation_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'FontSize', 9);
    
    % Add sample size information
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(results, gname)
            text(group_positions(g), mean_values(g) + sem_values(g) + 0.1*range(ylim), ...
                sprintf('n=%d', results.(gname).n_cells), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
    
    hold off;
    
    % Save figure
    saveas(gcf, fullfile(save_fig_folder, 'SpotSizeAnalysis_SI_BarChart.png'));
    saveas(gcf, fullfile(save_fig_folder, 'SpotSizeAnalysis_SI_BarChart.fig'));
    
    fprintf('\nFigure saved as: SpotSizeAnalysis_SI_BarChart.png\n');
else
    warning('No valid SI values calculated for any group.');
end

% ---------------- Individual Cell Traces Plot ----------------
% Create figure showing each cell trace of mean as function of spot size
if ~isempty(all_SI_values)
    figure('Color', 'w', 'Position', [200, 100, 1200, 600]);
    
    n_groups = numel(groups);
    
    for g = 1:n_groups
        gname = groups{g};
        if ~isfield(results, gname)
            continue;
        end
        
        % Create subplot for this group
        subplot(1, n_groups, g);
        
        % Get data for this group
        sizes = results.(gname).sizes;
        responses_matrix = results.(gname).responses_matrix;  % nSizes x nCells
        n_cells = size(responses_matrix, 2);
        
        % Plot each cell's response curve
        hold on;
        colors = lines(n_cells);  % Different color for each cell
        
        for c = 1:n_cells
            cell_responses = responses_matrix(:, c);
            
            % Only plot if cell has valid responses (not all NaN)
            if any(~isnan(cell_responses))
                plot(sizes, cell_responses, '-o', 'Color', colors(c,:), ...
                     'LineWidth', 1.5, 'MarkerSize', 6, 'MarkerFaceColor', colors(c,:), ...
                     'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            end
        end
        
        % Calculate and plot group mean
        mean_responses = mean(responses_matrix, 2, 'omitnan');
        sem_responses = std(responses_matrix, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(responses_matrix), 2));
        
        % Plot mean with error bars
        errorbar(sizes, mean_responses, sem_responses, 'k-', 'LineWidth', 3, ...
                 'MarkerSize', 8, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k');
        
        hold off;
        
        % Customize subplot
        xlabel('Spot Size (μm)');
        ylabel('Mean Response (spikes/s)');
        title(strrep(gname, '_', '\_'), 'FontSize', 12, 'FontWeight', 'bold');
        grid on;
        
        % Set x-axis to linear scale from 0 to 1200
        xlim([0, 1200]);
        set(gca, 'XScale', 'linear');
        xticks(0:200:1200);  % Tick marks every 200 μm
        
        % Add vertical lines to highlight size categories
        y_limits = ylim;
        
        % Mark large sizes used for S calculation
        for ls = large_sizes
            if any(sizes == ls)
                xline(ls, 'r--', 'LineWidth', 2, 'Alpha', 0.7);
            end
        end
        
        % Mark small size threshold
        xline(small_size_threshold, 'b--', 'LineWidth', 2, 'Alpha', 0.7);
        
        % Add legend for the first subplot
        if g == 1
            legend_entries = cell(n_cells + 1, 1);
            for c = 1:n_cells
                legend_entries{c} = sprintf('Cell %d', c);
            end
            legend_entries{end} = 'Group Mean ± SEM';
            legend(legend_entries, 'Location', 'best', 'FontSize', 8);
            
            % Add text box explaining line meanings
            text(0.02, 0.98, sprintf('Red dashed: Large sizes (%s)\nBlue dashed: Small size threshold (<%d)', ...
                mat2str(large_sizes), small_size_threshold), ...
                'Units', 'normalized', 'VerticalAlignment', 'top', ...
                'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'FontSize', 8);
        end
        
        % Add sample size to title
        title(sprintf('%s (n=%d)', strrep(gname, '_', '\_'), n_cells), ...
              'FontSize', 12, 'FontWeight', 'bold');
    end
    
    % Add overall title
    sgtitle(sprintf('Individual Cell Responses vs Spot Size (%s)', test_type), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, fullfile(save_fig_folder, sprintf('SpotSizeAnalysis_IndividualTraces_%s.png', test_type)));
    saveas(gcf, fullfile(save_fig_folder, sprintf('SpotSizeAnalysis_IndividualTraces_%s.fig', test_type)));
    
    fprintf('Individual traces figure saved as: SpotSizeAnalysis_IndividualTraces_%s.png\n', test_type);
end

% ---------------- Summary Statistics ----------------
fprintf('\n========== SIZE INDEX SUMMARY ==========\n');
for g = 1:numel(groups)
    gname = groups{g};
    if isfield(results, gname)
        fprintf('\n%s:\n', gname);
        fprintf('  n = %d cells\n', results.(gname).n_cells);
        fprintf('  Mean SI = %.3f ± %.3f (SEM)\n', results.(gname).mean_SI, results.(gname).sem_SI);
        fprintf('  Std SI = %.3f\n', results.(gname).std_SI);
        
        % Test if significantly different from 0
        valid_SI_values = results.(gname).SI_values(results.(gname).valid_SI);
        [~, p_val] = ttest(valid_SI_values, 0);
        fprintf('  t-test vs 0: p = %.4f', p_val);
        if p_val < 0.001
            fprintf(' ***\n');
        elseif p_val < 0.01
            fprintf(' **\n');
        elseif p_val < 0.05
            fprintf(' *\n');
        else
            fprintf(' (n.s.)\n');
        end
    end
end

% Statistical comparison between groups if we have exactly 2 groups
if numel(groups) == 2 && all(isfield(results, groups))
    fprintf('\n========== BETWEEN-GROUP COMPARISON ==========\n');
    group1_SI = results.(groups{1}).SI_values(results.(groups{1}).valid_SI);
    group2_SI = results.(groups{2}).SI_values(results.(groups{2}).valid_SI);
    
    [~, p_val_between] = ttest2(group1_SI, group2_SI);
    fprintf('t-test between %s and %s: p = %.4f', groups{1}, groups{2}, p_val_between);
    if p_val_between < 0.001
        fprintf(' ***\n');
    elseif p_val_between < 0.01
        fprintf(' **\n');
    elseif p_val_between < 0.05
        fprintf(' *\n');
    else
        fprintf(' (n.s.)\n');
    end
end

% ---------------- Optimal Spot Size Analysis ----------------
% Calculate optimal spot size (specified percentage of peak firing rate) for each group
fprintf('\n========== OPTIMAL SPOT SIZE ANALYSIS (%.0f%% of Peak) ==========\n', percent_peak*100);

optimal_results = struct();
all_optimal_sizes = [];
all_optimal_group_labels = [];

for g = 1:numel(groups)
    gname = groups{g};
    if ~isfield(results, gname)
        continue;
    end
    
    fprintf('\nProcessing optimal spot size for group %s:\n', gname);
    
    % Get data for this group
    sizes = results.(gname).sizes;
    responses_matrix = results.(gname).responses_matrix;  % nSizes x nCells
    n_cells = size(responses_matrix, 2);
    
    optimal_sizes = nan(n_cells, 1);
    
    for c = 1:n_cells
        cell_responses = responses_matrix(:, c);
        
        % Skip if cell has too many NaN values
        valid_responses = ~isnan(cell_responses);
        if sum(valid_responses) < 3
            continue;
        end
        
        % Find peak response
        peak_response = max(cell_responses, [], 'omitnan');
        
        % Calculate specified percentage of peak
        target_response = percent_peak * peak_response;
        
        % Find the spot size that first reaches the target percentage of peak (ascending order)
        target_reached = cell_responses >= target_response;
        first_target_idx = find(target_reached, 1, 'first');
        
        if ~isempty(first_target_idx)
            optimal_sizes(c) = sizes(first_target_idx);
            
            % Debug output for first few cells
            if c <= 3
                fprintf('  Cell %d: Peak=%.2f, %.0f%% target=%.2f, Optimal size=%d\n', ...
                    c, peak_response, percent_peak*100, target_response, optimal_sizes(c));
            end
        end
    end
    
    % Remove NaN values
    valid_optimal = ~isnan(optimal_sizes);
    valid_optimal_sizes = optimal_sizes(valid_optimal);
    
    fprintf('  Valid optimal size calculations: %d/%d cells\n', sum(valid_optimal), n_cells);
    if ~isempty(valid_optimal_sizes)
        fprintf('  Optimal size range: %d to %d μm\n', min(valid_optimal_sizes), max(valid_optimal_sizes));
        fprintf('  Mean optimal size: %.1f ± %.1f μm (SEM)\n', mean(valid_optimal_sizes), std(valid_optimal_sizes)/sqrt(length(valid_optimal_sizes)));
        
        % Store results
        optimal_results.(gname).optimal_sizes = optimal_sizes;
        optimal_results.(gname).valid_optimal = valid_optimal;
        optimal_results.(gname).mean_optimal = mean(valid_optimal_sizes);
        optimal_results.(gname).sem_optimal = std(valid_optimal_sizes) / sqrt(length(valid_optimal_sizes));
        optimal_results.(gname).std_optimal = std(valid_optimal_sizes);
        optimal_results.(gname).n_cells = length(valid_optimal_sizes);
        
        % Collect data for combined plot
        all_optimal_sizes = [all_optimal_sizes; valid_optimal_sizes];
        all_optimal_group_labels = [all_optimal_group_labels; repmat(g, length(valid_optimal_sizes), 1)];
    else
        fprintf('  No valid optimal sizes calculated\n');
    end
end

% ---------------- Optimal Spot Size Plotting ----------------
if ~isempty(all_optimal_sizes)
    figure('Color', 'w', 'Position', [300, 100, 800, 600]);
    
    % Calculate bar positions
    group_positions = 1:numel(groups);
    bar_width = 0.6;
    
    % Prepare data for bar plot
    mean_optimal_values = nan(numel(groups), 1);
    sem_optimal_values = nan(numel(groups), 1);
    
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(optimal_results, gname)
            mean_optimal_values(g) = optimal_results.(gname).mean_optimal;
            sem_optimal_values(g) = optimal_results.(gname).sem_optimal;
        end
    end
    
    % Create bar plot with error bars
    hold on;
    bars = bar(group_positions, mean_optimal_values, bar_width, 'FaceAlpha', 0.7, 'FaceColor', [0.5 0.8 0.5]);
    
    % Add error bars
    errorbar(group_positions, mean_optimal_values, sem_optimal_values, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'CapSize', 8);
    
    % Add individual data points
    jitter_width = 0.15;
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(optimal_results, gname)
            valid_optimal_values = optimal_results.(gname).optimal_sizes(optimal_results.(gname).valid_optimal);
            
            % Add jitter to x-positions
            n_points = length(valid_optimal_values);
            x_jitter = (rand(n_points, 1) - 0.5) * jitter_width * 2;
            x_positions = repmat(group_positions(g), n_points, 1) + x_jitter;
            
            % Plot individual points
            scatter(x_positions, valid_optimal_values, 40, [0.2 0.6 0.2], 'filled', ...
                'MarkerEdgeColor', 'k', 'LineWidth', 0.5, 'MarkerFaceAlpha', 0.8);
        end
    end
    
    % Customize plot
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    xlabel('Cell Groups');
    ylabel('Optimal Spot Size (μm)');
    title(sprintf('Optimal Spot Size (%.0f%% of Peak Response)', percent_peak*100));
    grid on;
    
    % Add text annotation explaining calculation
    annotation_text = sprintf('Optimal size = smallest spot size\nthat reaches %.0f%% of peak firing rate', percent_peak*100);
    text(0.02, 0.98, annotation_text, 'Units', 'normalized', ...
        'VerticalAlignment', 'top', 'HorizontalAlignment', 'left', ...
        'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k', 'FontSize', 9);
    
    % Add sample size information
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(optimal_results, gname)
            text(group_positions(g), mean_optimal_values(g) + sem_optimal_values(g) + 0.05*range(ylim), ...
                sprintf('n=%d', optimal_results.(gname).n_cells), ...
                'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
        end
    end
    
    hold off;
    
    % Save figure
    saveas(gcf, fullfile(save_fig_folder, sprintf('OptimalSpotSize_%s.png', test_type)));
    saveas(gcf, fullfile(save_fig_folder, sprintf('OptimalSpotSize_%s.fig', test_type)));
    
    fprintf('\nOptimal spot size figure saved as: OptimalSpotSize_%s.png\n', test_type);
end

% ---------------- Optimal Spot Size Statistics ----------------
fprintf('\n========== OPTIMAL SPOT SIZE SUMMARY (%.0f%% of Peak) ==========\n', percent_peak*100);
for g = 1:numel(groups)
    gname = groups{g};
    if isfield(optimal_results, gname)
        fprintf('\n%s:\n', gname);
        fprintf('  n = %d cells\n', optimal_results.(gname).n_cells);
        fprintf('  Mean optimal size = %.1f ± %.1f μm (SEM)\n', optimal_results.(gname).mean_optimal, optimal_results.(gname).sem_optimal);
        fprintf('  Std optimal size = %.1f μm\n', optimal_results.(gname).std_optimal);
        
        % Test distribution of optimal sizes
        valid_optimal_values = optimal_results.(gname).optimal_sizes(optimal_results.(gname).valid_optimal);
        fprintf('  Range: %d - %d μm\n', min(valid_optimal_values), max(valid_optimal_values));
        fprintf('  Median: %.1f μm\n', median(valid_optimal_values));
    end
end

% Statistical comparison between groups for optimal sizes
if numel(groups) == 2 && all(isfield(optimal_results, groups))
    fprintf('\n========== OPTIMAL SIZE BETWEEN-GROUP COMPARISON ==========\n');
    group1_optimal = optimal_results.(groups{1}).optimal_sizes(optimal_results.(groups{1}).valid_optimal);
    group2_optimal = optimal_results.(groups{2}).optimal_sizes(optimal_results.(groups{2}).valid_optimal);
    
    [~, p_val_optimal] = ttest2(group1_optimal, group2_optimal);
    fprintf('t-test between %s and %s: p = %.4f', groups{1}, groups{2}, p_val_optimal);
    if p_val_optimal < 0.001
        fprintf(' ***\n');
    elseif p_val_optimal < 0.01
        fprintf(' **\n');
    elseif p_val_optimal < 0.05
        fprintf(' *\n');
    else
        fprintf(' (n.s.)\n');
    end
    
    % Also perform Mann-Whitney U test (non-parametric)
    p_val_mw = ranksum(group1_optimal, group2_optimal);
    fprintf('Mann-Whitney U test: p = %.4f', p_val_mw);
    if p_val_mw < 0.001
        fprintf(' ***\n');
    elseif p_val_mw < 0.01
        fprintf(' **\n');
    elseif p_val_mw < 0.05
        fprintf(' *\n');
    else
        fprintf(' (n.s.)\n');
    end
end

fprintf('\nAnalysis complete!\n');
