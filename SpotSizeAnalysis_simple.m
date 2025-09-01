% Simplified spot size analysis: Calculate Size Index (SI) for different groups
% SI = (S-C)/C where:
%   S = mean response of spot sizes 800, 1200
%   C = mean response of top two spot sizes < 800
%
% Requires: struct 'a' in workspace with fields:
%   - x1 (T x 1), y_labels (1 x S), and group arrays (T x (S * Ncells))
% Example groups: 'DN_ONSus_RF_UV', 'DN_ONSus_RF_GRN'

clc; close all;
test_type = 'ON';
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
grey_color = 0.3*ones(1, 3);

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
% Plot id: 12943
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
    bars = bar(group_positions, mean_values, bar_width, 'FaceAlpha', 0.7, 'EdgeColor', 'none');
    
    % Add error bars
    errorbar(group_positions, mean_values, sem_values, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'CapSize', 0);
    
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
            scatter(x_positions, valid_SI_values, 40, grey_color, 'filled', ...
                'MarkerEdgeColor', 'none');
        end
    end
    
    % Customize plot
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    ylim([-0.25 0.1])
    yticks(-0.2:0.1:0.1)
    yticklabels({'-0.2', '-0.1', '0', '0.1'})
    xlabel('Cell Groups');
    ylabel('Size Index (SI)');
    title('Size Index: (Large Spot Response - Center Response) / Center Response');
    grid off;
    
    % Add horizontal line at SI = 0
    yline(0, 'k--', 'LineWidth', 1);
    
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
    print(gcf, fullfile(save_fig_folder, 'SpotSizeAnalysis_SI_BarChart.eps'), '-depsc', '-r300');
    
    fprintf('\nFigure saved as: SpotSizeAnalysis_SI_BarChart.png and .eps (vector)\n');
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
                xline(ls, 'r--', 'LineWidth', 2);
            end
        end
        
        % Mark small size threshold
        xline(small_size_threshold, 'b--', 'LineWidth', 2);
        
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
    
    fprintf('Individual traces figure saved as: SpotSizeAnalysis_IndividualTraces_%s.png\n', test_type);
end

% ---------------- C and S Temporal Traces Plot ----------------
% Create figure showing temporal traces for C and S spot sizes
if ~isempty(all_SI_values)
    fprintf('\n========== CREATING C AND S TEMPORAL TRACES ==========\n');
    
    figure('Color', 'w', 'Position', [400, 100, 1200, 800]);
    
    n_groups = numel(groups);
    
    for g = 1:n_groups
        gname = groups{g};
        if ~isfield(results, gname)
            continue;
        end
        
        % Create subplot for this group
        subplot(2, n_groups, g);
        
        % Get original temporal data for this group
        M = a.(gname);            % T x (nSizes * nCells)
        [Tg, nCol] = size(M);
        nCells = nCol / nSizes;
        
        % Reshape temporal data: T x nSizes x nCells
        M_reshaped = reshape(M, [Tg, nSizes, nCells]);
        
        % Reorder sizes to be ascending (and reorder M_reshaped accordingly)
        M_reshaped = M_reshaped(:, sort_idx, :);
        
        % Initialize arrays to store C and S spot indices for each cell
        C_spot_indices = cell(nCells, 1);
        S_spot_indices = cell(nCells, 1);
        
        % Re-calculate C and S indices for each cell (following exact same logic as SI calculation)
        for c = 1:nCells
            cell_responses = squeeze(mean(M_reshaped(stim_idx, :, c), 1, 'omitnan'));  % Mean during stim for this cell
            
            % Find indices for large sizes (S) - same as in SI calculation
            large_size_indices = ismember(sorted_sizes, large_sizes);
            if sum(large_size_indices) > 0
                S_spot_indices{c} = find(large_size_indices);
            end
            
            % Find indices for small sizes (C) - same as in SI calculation  
            small_size_indices = sorted_sizes < small_size_threshold;
            if sum(small_size_indices) >= 2
                C_responses = cell_responses(small_size_indices);
                small_sizes_for_this_cell = sorted_sizes(small_size_indices);
                
                % Remove NaN values for peak finding
                valid_indices = ~isnan(C_responses);
                if sum(valid_indices) >= 2
                    valid_C_responses = C_responses(valid_indices);
                    valid_small_sizes = small_sizes_for_this_cell(valid_indices);
                    
                    % Find peak response
                    [~, peak_idx] = max(valid_C_responses);
                    peak_size = valid_small_sizes(peak_idx);
                    
                    % Find peak size index in the original sorted_sizes array
                    peak_size_original_idx = find(sorted_sizes == peak_size);
                    
                    % Find adjacent spot sizes to the peak (same logic as SI calculation)
                    adjacent_responses = [];
                    adjacent_indices = [];
                    
                    % Check size immediately before peak
                    if peak_size_original_idx > 1
                        prev_size = sorted_sizes(peak_size_original_idx - 1);
                        if prev_size < small_size_threshold && ~isnan(cell_responses(peak_size_original_idx - 1))
                            adjacent_responses(end+1) = cell_responses(peak_size_original_idx - 1);
                            adjacent_indices(end+1) = peak_size_original_idx - 1;
                        end
                    end
                    
                    % Check size immediately after peak
                    if peak_size_original_idx < length(sorted_sizes)
                        next_size = sorted_sizes(peak_size_original_idx + 1);
                        if next_size < small_size_threshold && ~isnan(cell_responses(peak_size_original_idx + 1))
                            adjacent_responses(end+1) = cell_responses(peak_size_original_idx + 1);
                            adjacent_indices(end+1) = peak_size_original_idx + 1;
                        end
                    end
                    
                    % Store C indices: peak + largest adjacent (same as SI calculation)
                    if length(adjacent_responses) >= 1
                        [~, largest_idx] = max(adjacent_responses);
                        C_spot_indices{c} = [peak_size_original_idx, adjacent_indices(largest_idx)];
                        
                        % Debug for first few cells
                        if c <= 3
                            fprintf('  Group %s, Cell %d: C spots = %s (sizes: %s), S spots = %s (sizes: %s)\n', ...
                                gname, c, mat2str(C_spot_indices{c}), mat2str(sorted_sizes(C_spot_indices{c})), ...
                                mat2str(S_spot_indices{c}), mat2str(sorted_sizes(S_spot_indices{c})));
                        end
                    end
                end
            end
        end
        
        % Calculate mean temporal traces for C and S across all cells
        C_traces = [];
        S_traces = [];
        
        for c = 1:nCells
            % Get C trace for this cell (mean across C spot sizes)
            if ~isempty(C_spot_indices{c})
                C_trace_cell = squeeze(mean(M_reshaped(:, C_spot_indices{c}, c), 2, 'omitnan'));
                C_traces = [C_traces, C_trace_cell];
            end
            
            % Get S trace for this cell (mean across S spot sizes)
            if ~isempty(S_spot_indices{c})
                S_trace_cell = squeeze(mean(M_reshaped(:, S_spot_indices{c}, c), 2, 'omitnan'));
                S_traces = [S_traces, S_trace_cell];
            end
        end
        
        % Calculate group means and SEMs
        if ~isempty(C_traces)
            C_mean = mean(C_traces, 2, 'omitnan');
            C_sem = std(C_traces, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(C_traces), 2));
        else
            C_mean = nan(Tg, 1);
            C_sem = nan(Tg, 1);
        end
        
        if ~isempty(S_traces)
            S_mean = mean(S_traces, 2, 'omitnan');
            S_sem = std(S_traces, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(S_traces), 2));
        else
            S_mean = nan(Tg, 1);
            S_sem = nan(Tg, 1);
        end
        
        % Plot temporal traces
        hold on;
        time_axis = a.x1;
        
        % Plot C trace (center) in blue
        if ~isempty(C_traces)
            fill([time_axis; flipud(time_axis)], [C_mean + C_sem; flipud(C_mean - C_sem)], ...
                 [0.3 0.3 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(time_axis, C_mean, 'b-', 'LineWidth', 2, 'DisplayName', 'C (Center)');
        end
        
        % Plot S trace (surround) in red
        if ~isempty(S_traces)
            fill([time_axis; flipud(time_axis)], [S_mean + S_sem; flipud(S_mean - S_sem)], ...
                 [1 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(time_axis, S_mean, 'r-', 'LineWidth', 2, 'DisplayName', 'S (Surround)');
        end
        
        % Add stimulus period shading
        y_lim = ylim;
        if strcmp(test_type, 'ON')
            fill([time_axis(110), time_axis(300), time_axis(300), time_axis(110)], ...
                 [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        elseif strcmp(test_type, 'OFF')
            fill([time_axis(310), time_axis(500), time_axis(500), time_axis(310)], ...
                 [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], [0.8 0.8 0.8], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end
        
        hold off;
        
        % Customize subplot
        xlabel('Time (s)');
        ylabel('Response (spikes/s)');
        title(sprintf('%s\nC vs S Temporal Traces', strrep(gname, '_', '\_')), ...
              'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        legend('Location', 'best');
        
        % Add sample sizes to title
        n_C_cells = sum(~cellfun(@isempty, C_spot_indices));
        n_S_cells = sum(~cellfun(@isempty, S_spot_indices));
        title(sprintf('%s\nC vs S Temporal Traces (C: n=%d, S: n=%d)', ...
              strrep(gname, '_', '\_'), n_C_cells, n_S_cells), ...
              'FontSize', 11, 'FontWeight', 'bold');
        
        % Second subplot: Show stimulus period only (zoomed in)
        subplot(2, n_groups, g + n_groups);
        
        % Define stimulus period indices
        if strcmp(test_type, 'ON')
            stim_period_idx = 110:300;
        elseif strcmp(test_type, 'OFF')
            stim_period_idx = 310:500;
        end
        
        hold on;
        time_axis_stim = time_axis(stim_period_idx);
        
        % Plot C trace during stimulus period
        if ~isempty(C_traces)
            C_mean_stim = C_mean(stim_period_idx);
            C_sem_stim = C_sem(stim_period_idx);
            fill([time_axis_stim; flipud(time_axis_stim)], [C_mean_stim + C_sem_stim; flipud(C_mean_stim - C_sem_stim)], ...
                 [0.3 0.3 1], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(time_axis_stim, C_mean_stim, 'b-', 'LineWidth', 2, 'DisplayName', 'C (Center)');
        end
        
        % Plot S trace during stimulus period
        if ~isempty(S_traces)
            S_mean_stim = S_mean(stim_period_idx);
            S_sem_stim = S_sem(stim_period_idx);
            fill([time_axis_stim; flipud(time_axis_stim)], [S_mean_stim + S_sem_stim; flipud(S_mean_stim - S_sem_stim)], ...
                 [1 0.3 0.3], 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            plot(time_axis_stim, S_mean_stim, 'r-', 'LineWidth', 2, 'DisplayName', 'S (Surround)');
        end
        
        hold off;
        
        % Customize stimulus period subplot
        xlabel('Time (s)');
        ylabel('Response (spikes/s)');
        title(sprintf('Stimulus Period (%s)', test_type), 'FontSize', 11);
        grid on;
        legend('Location', 'best');
        
        fprintf('  Group %s: C traces from %d cells, S traces from %d cells\n', gname, n_C_cells, n_S_cells);
    end
    
    % Add overall title
    sgtitle(sprintf('Temporal Traces for C (Center) and S (Surround) Spot Sizes (%s)', test_type), ...
            'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, fullfile(save_fig_folder, sprintf('C_S_TemporalTraces_%s.png', test_type)));
    print(gcf, fullfile(save_fig_folder, sprintf('C_S_TemporalTraces_%s.eps', test_type)), '-depsc', '-r300');
    
    fprintf('\nC and S temporal traces figure saved as: C_S_TemporalTraces_%s.png and .eps\n', test_type);
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

% Parameters for interpolation
interp_resolution = 10; % μm spacing for interpolation
interp_sizes = min(sorted_sizes):interp_resolution:max(sorted_sizes);

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
    interpolated_curves = cell(n_cells, 1); % Store interpolated curves for plotting
    
    for c = 1:n_cells
        cell_responses = responses_matrix(:, c);
        
        % Skip if cell has too many NaN values
        valid_idx = ~isnan(cell_responses);
        if sum(valid_idx) < 3
            continue;
        end
        
        valid_sizes = sizes(valid_idx);
        valid_responses = cell_responses(valid_idx);
        
        % Use pchip interpolation (shape-preserving, no overshoot)
        try
            % Create interpolation function
            pp = pchip(valid_sizes, valid_responses);
            interp_responses = ppval(pp, interp_sizes);
            
            % Store interpolated curve
            interpolated_curves{c} = struct('sizes', interp_sizes, 'responses', interp_responses);
            
            % Find peak response from interpolated data
            peak_response = max(interp_responses);
            
            % Calculate specified percentage of peak
            target_response = percent_peak * peak_response;
            
            % Find the spot size that first reaches the target percentage of peak
            target_reached = interp_responses >= target_response;
            first_target_idx = find(target_reached, 1, 'first');
            
            if ~isempty(first_target_idx)
                optimal_sizes(c) = interp_sizes(first_target_idx);
                
                % Debug output for first few cells
                if c <= 3
                    fprintf('  Cell %d: Peak=%.2f, %.0f%% target=%.2f, Optimal size=%.1f μm\n', ...
                        c, peak_response, percent_peak*100, target_response, optimal_sizes(c));
                end
            end
        catch ME
            fprintf('  Cell %d: Interpolation failed - %s\n', c, ME.message);
        end
    end
    
    % Remove NaN values
    valid_optimal = ~isnan(optimal_sizes);
    valid_optimal_sizes = optimal_sizes(valid_optimal);
    
    fprintf('  Valid optimal size calculations: %d/%d cells\n', sum(valid_optimal), n_cells);
    if ~isempty(valid_optimal_sizes)
        fprintf('  Optimal size range: %.1f to %.1f μm\n', min(valid_optimal_sizes), max(valid_optimal_sizes));
        fprintf('  Mean optimal size: %.1f ± %.1f μm (SEM)\n', mean(valid_optimal_sizes), std(valid_optimal_sizes)/sqrt(length(valid_optimal_sizes)));
        
        % Store results
        optimal_results.(gname).optimal_sizes = optimal_sizes;
        optimal_results.(gname).valid_optimal = valid_optimal;
        optimal_results.(gname).mean_optimal = mean(valid_optimal_sizes);
        optimal_results.(gname).sem_optimal = std(valid_optimal_sizes) / sqrt(length(valid_optimal_sizes));
        optimal_results.(gname).std_optimal = std(valid_optimal_sizes);
        optimal_results.(gname).n_cells = length(valid_optimal_sizes);
        optimal_results.(gname).interpolated_curves = interpolated_curves;
        
        % Collect data for combined plot
        all_optimal_sizes = [all_optimal_sizes; valid_optimal_sizes];
        all_optimal_group_labels = [all_optimal_group_labels; repmat(g, length(valid_optimal_sizes), 1)];
    else
        fprintf('  No valid optimal sizes calculated\n');
    end
end

% ---------------- Interpolation Visualization ----------------
% Create figure showing interpolated curves for each group
if ~isempty(all_optimal_sizes)
    figure('Color', 'w', 'Position', [100, 50, 1400, 800]);
    
    n_groups = numel(groups);
    
    for g = 1:n_groups
        gname = groups{g};
        if ~isfield(optimal_results, gname)
            continue;
        end
        
        % Create subplot for this group
        subplot(2, n_groups, g);
        
        % Get data for this group
        sizes = results.(gname).sizes;
        responses_matrix = results.(gname).responses_matrix;  % nSizes x nCells
        interpolated_curves = optimal_results.(gname).interpolated_curves;
        optimal_sizes_for_group = optimal_results.(gname).optimal_sizes;
        n_cells = size(responses_matrix, 2);
        
        % Plot each cell's interpolated curve
        hold on;
        colors = lines(n_cells);
        
        for c = 1:n_cells
            cell_responses = responses_matrix(:, c);
            
            % Plot original data points
            valid_idx = ~isnan(cell_responses);
            if any(valid_idx)
                plot(sizes(valid_idx), cell_responses(valid_idx), 'o', ...
                     'Color', colors(c,:), 'MarkerSize', 6, ...
                     'MarkerFaceColor', colors(c,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            end
            
            % Plot interpolated curve if available
            if ~isempty(interpolated_curves{c})
                plot(interpolated_curves{c}.sizes, interpolated_curves{c}.responses, '-', ...
                     'Color', colors(c,:), 'LineWidth', 1.5);
                
                % Mark optimal spot size if valid
                if ~isnan(optimal_sizes_for_group(c))
                    % Find corresponding response at optimal size
                    [~, opt_idx] = min(abs(interpolated_curves{c}.sizes - optimal_sizes_for_group(c)));
                    opt_response = interpolated_curves{c}.responses(opt_idx);
                    
                    plot(optimal_sizes_for_group(c), opt_response, 's', ...
                         'Color', colors(c,:), 'MarkerSize', 10, ...
                         'MarkerFaceColor', 'none', 'MarkerEdgeColor', colors(c,:), 'LineWidth', 2);
                end
            end
        end
        
        hold off;
        
        % Customize subplot
        xlabel('Spot Size (μm)');
        ylabel('Mean Response (spikes/s)');
        title(sprintf('%s\nInterpolated Curves (pchip)', strrep(gname, '_', '\_')), ...
              'FontSize', 11, 'FontWeight', 'bold');
        grid on;
        xlim([0, 1200]);
        set(gca, 'XScale', 'linear');
        xticks(0:200:1200);
        
        % Add legend for the first subplot
        if g == 1
            legend_str = sprintf('Circles: Original data\nLines: Interpolated\nSquares: %.0f%% optimal', percent_peak*100);
            text(0.02, 0.98, legend_str, 'Units', 'normalized', ...
                 'VerticalAlignment', 'top', 'BackgroundColor', [1 1 1 0.8], ...
                 'EdgeColor', 'k', 'FontSize', 8);
        end
        
        % Lower subplot: Show zoom-in around optimal sizes
        subplot(2, n_groups, g + n_groups);
        
        hold on;
        
        % Determine zoom range based on optimal sizes
        valid_optimal = optimal_sizes_for_group(~isnan(optimal_sizes_for_group));
        if ~isempty(valid_optimal)
            zoom_center = mean(valid_optimal);
            zoom_range = max(200, 2 * std(valid_optimal));
            zoom_min = max(0, zoom_center - zoom_range);
            zoom_max = min(1200, zoom_center + zoom_range);
        else
            zoom_min = 0;
            zoom_max = 600;
        end
        
        for c = 1:n_cells
            cell_responses = responses_matrix(:, c);
            
            % Plot original data points in zoom range
            valid_idx = ~isnan(cell_responses) & sizes >= zoom_min & sizes <= zoom_max;
            if any(valid_idx)
                plot(sizes(valid_idx), cell_responses(valid_idx), 'o', ...
                     'Color', colors(c,:), 'MarkerSize', 8, ...
                     'MarkerFaceColor', colors(c,:), 'MarkerEdgeColor', 'k', 'LineWidth', 0.5);
            end
            
            % Plot interpolated curve in zoom range
            if ~isempty(interpolated_curves{c})
                zoom_idx = interpolated_curves{c}.sizes >= zoom_min & ...
                          interpolated_curves{c}.sizes <= zoom_max;
                plot(interpolated_curves{c}.sizes(zoom_idx), ...
                     interpolated_curves{c}.responses(zoom_idx), '-', ...
                     'Color', colors(c,:), 'LineWidth', 2);
                
                % Mark optimal spot size
                if ~isnan(optimal_sizes_for_group(c))
                    % Find peak response for this cell
                    peak_response = max(interpolated_curves{c}.responses);
                    target_response = percent_peak * peak_response;
                    
                    % Plot horizontal line at target response
                    plot([zoom_min, zoom_max], [target_response, target_response], '--', ...
                         'Color', colors(c,:), 'LineWidth', 1);
                    
                    % Mark optimal point
                    [~, opt_idx] = min(abs(interpolated_curves{c}.sizes - optimal_sizes_for_group(c)));
                    opt_response = interpolated_curves{c}.responses(opt_idx);
                    
                    plot(optimal_sizes_for_group(c), opt_response, 's', ...
                         'Color', colors(c,:), 'MarkerSize', 12, ...
                         'MarkerFaceColor', colors(c,:), 'MarkerEdgeColor', 'k', 'LineWidth', 2);
                    
                    % Add vertical line at optimal size
                    plot([optimal_sizes_for_group(c), optimal_sizes_for_group(c)], ...
                         [0, opt_response], ':', 'Color', colors(c,:), 'LineWidth', 1.5);
                end
            end
        end
        
        % Add mean optimal size line
        if ~isempty(valid_optimal)
            xline(mean(valid_optimal), 'k-', 'LineWidth', 2, ...
                  'Label', sprintf('Mean: %.1f μm', mean(valid_optimal)));
        end
        
        hold off;
        
        % Customize zoom subplot
        xlabel('Spot Size (μm)');
        ylabel('Mean Response (spikes/s)');
        title(sprintf('Zoomed View (n=%d cells)', sum(~isnan(optimal_sizes_for_group))), ...
              'FontSize', 11);
        grid on;
        xlim([zoom_min, zoom_max]);
        set(gca, 'XScale', 'linear');
    end
    
    % Add overall title
    sgtitle(sprintf('Interpolated Response Curves and Optimal Spot Sizes (%.0f%% of Peak) - %s', ...
            percent_peak*100, test_type), 'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, fullfile(save_fig_folder, sprintf('OptimalSpotSize_Interpolated_%s.png', test_type)));
    
    fprintf('\nInterpolated curves figure saved as: OptimalSpotSize_Interpolated_%s.png\n', test_type);
end

% ---------------- Optimal Spot Size Plotting ----------------
% plot_id: 48521
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
    bars = bar(group_positions, mean_optimal_values, bar_width, 'FaceAlpha', 0.7, 'FaceColor', [0.5 0.8 0.5], 'EdgeColor', 'none');
    
    % Add error bars
    errorbar(group_positions, mean_optimal_values, sem_optimal_values, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'CapSize', 0);
    
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
            scatter(x_positions, valid_optimal_values, 40, grey_color, 'filled', ...
                'MarkerEdgeColor', 'none');
        end
    end
    
    % Customize plot
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    ylim([0, 250]);
    yticks(0:100:200);
    yticklabels({'0', '100', '200'});
    xlabel('Cell Groups');
    ylabel('Optimal Spot Size (μm)');
    title(sprintf('Optimal Spot Size (%.0f%% of Peak Response) - Interpolated', percent_peak*100));
    grid off;

    % Add text annotation explaining calculation
    annotation_text = sprintf(['Optimal size = smallest spot size\n' ...
                              'that reaches %.0f%% of peak firing rate\n' ...
                              'Using pchip interpolation (%.0f μm resolution)'], ...
                              percent_peak*100, interp_resolution);
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
    saveas(gcf, fullfile(save_fig_folder, sprintf('OptimalSpotSize_Interpolated_BarChart_%s.png', test_type)));
    print(gcf, fullfile(save_fig_folder, sprintf('OptimalSpotSize_Interpolated_BarChart_%s.eps', test_type)), '-depsc', '-r300');
    
    fprintf('\nOptimal spot size bar chart saved as: OptimalSpotSize_Interpolated_BarChart_%s.png and .eps (vector)\n', test_type);
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
