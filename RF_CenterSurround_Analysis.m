% RF_CenterSurround_Analysis.m
% Estimate center-surround receptive field parameters using estimate_rf_sizes function
% Analyzes both individual RGCs and group averages
%
% Requires:
%   - estimate_rf_sizes.m in path
%   - Data struct 'a' with spot size response data
%
% Author: Created for PreyCaptureRGC analysis
% Date: September 10, 2025

clc; close all;

% ================== CONFIGURATION ==================
test_type = 'ON';  % 'ON' or 'OFF'
draws = 2000;      % Number of bootstrap draws for RF estimation
bigR = [400, 600, 800];  % Large radii for surround estimation (μm)

% Save results folder
save_folder = './Results/RF_Analysis/';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end

% Load data
fprintf('Loading data...\n');
if ~exist('a', 'var')
    a = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\Spots\SusAlpha_RawData.mat');
end

% ================== SETUP GROUPS AND PARAMETERS ==================
switch test_type
    case 'ON'
        groups = {'AcuteZoneDT_ONSus_RF_GRN', 'DN_ONSus_RF_GRN', ...
                  'AcuteZoneDT_ONSus_RF_UV', 'DN_ONSus_RF_UV'};
        stim_idx = 110:300;  % during-stimulus period for ON responses
    case 'OFF'  
        groups = {'AcuteZoneDT_OFFSus_RF_GRN', 'DN_OFFSus_RF_GRN', ...
                  'AcuteZoneDT_OFFSus_RF_UV', 'DN_OFFSus_RF_UV'};
        stim_idx = 310:500;  % during-stimulus period for OFF responses
end

% Get spot sizes (convert from diameters to radii)
sizes_diameter = a.y_labels(:)';  % spot diameters in μm
sizes_radius = sizes_diameter / 2;  % convert to radii for RF estimation
nSizes = numel(sizes_diameter);
T = numel(a.x1);

% Sort sizes
[sorted_radii, sort_idx] = sort(sizes_radius, 'ascend');
sorted_diameters = sorted_radii * 2;

fprintf('Spot sizes (diameters): %s μm\n', mat2str(sorted_diameters));
fprintf('Spot sizes (radii): %s μm\n', mat2str(sorted_radii));
fprintf('Bootstrap draws: %d\n', draws);
fprintf('Large radii for surround estimation: %s μm\n', mat2str(bigR));

% ================== COLOR CONFIGURATION ==================
group_colors = containers.Map();
group_colors('AcuteZoneDT_ONSus_RF_GRN') = [180, 0, 180]/255;   % Magenta
group_colors('DN_ONSus_RF_GRN') = [120, 0, 120]/255;           % Dark Magenta  
group_colors('AcuteZoneDT_OFFSus_RF_GRN') = [0, 180, 0]/255;   % Green
group_colors('DN_OFFSus_RF_GRN') = [0, 120, 0]/255;            % Dark Green
group_colors('AcuteZoneDT_ONSus_RF_UV') = [180, 180, 180]/255; % Light Gray
group_colors('DN_ONSus_RF_UV') = [120, 120, 120]/255;          % Dark Gray
group_colors('AcuteZoneDT_OFFSus_RF_UV') = [180, 180, 180]/255; % Light Gray
group_colors('DN_OFFSus_RF_UV') = [120, 120, 120]/255;         % Dark Gray

% ================== MAIN ANALYSIS ==================
rf_results = struct();
all_individual_results = [];
all_group_labels = [];

fprintf('\n========== RF CENTER-SURROUND ANALYSIS ==========\n');

for g = 1:numel(groups)
    gname = groups{g};
    
    if ~isfield(a, gname)
        warning('Group "%s" not found in struct a. Skipping.', gname);
        continue;
    end
    
    fprintf('\n--- Processing group: %s ---\n', gname);
    
    % Get data for this group
    M = a.(gname);  % T x (nSizes * nCells)
    [Tg, nCol] = size(M);
    
    % Validate dimensions
    if Tg ~= T
        error('Time length mismatch for group %s', gname);
    end
    if mod(nCol, nSizes) ~= 0
        error('Columns not divisible by number of sizes for group %s', gname);
    end
    
    nCells = nCol / nSizes;
    fprintf('  Number of cells: %d\n', nCells);
    
    % Calculate mean response during stimulus for each cell-size combination
    stim_mean_by_col = mean(M(stim_idx, :), 1, 'omitnan');
    
    % Reshape to [nSizes x nCells] and sort by size
    response_matrix = reshape(stim_mean_by_col, [nSizes, nCells]);
    response_matrix = response_matrix(sort_idx, :);  % Sort by ascending size
    
    % ============== INDIVIDUAL CELL RF ANALYSIS ==============
    fprintf('  Analyzing individual cells...\n');
    
    individual_results = struct();
    individual_results.group = gname;
    individual_results.n_cells = nCells;
    individual_results.rf_params = cell(nCells, 1);
    individual_results.valid_cells = false(nCells, 1);
    
    % Storage for individual cell parameters
    sigma_c_ind = nan(nCells, 1);
    sigma_s_ind = nan(nCells, 1);
    diam63_c_ind = nan(nCells, 1);
    diam63_s_ind = nan(nCells, 1);
    
    for c = 1:nCells
        cell_responses = response_matrix(:, c);
        
        % Skip cells with too many NaN values or insufficient data
        valid_responses = ~isnan(cell_responses) & cell_responses >= 0;
        if sum(valid_responses) < 4  % Need at least 4 points for fitting
            fprintf('    Cell %d: Insufficient valid data points\n', c);
            continue;
        end
        
        % Prepare data for RF estimation
        rvals = sorted_radii(valid_responses);
        Rmean = cell_responses(valid_responses)';
        Rsem = [];  % No SEM for individual cells
        
        % Skip if no clear peak or all responses are too similar
        if max(Rmean) - min(Rmean) < 0.1 * max(Rmean)
            fprintf('    Cell %d: Insufficient response modulation\n', c);
            continue;
        end
        
        try
            % Estimate RF parameters
            rf_params = estimate_rf_sizes(rvals, Rmean, Rsem, bigR, draws);
            
            % Store results
            individual_results.rf_params{c} = rf_params;
            individual_results.valid_cells(c) = true;
            
            sigma_c_ind(c) = rf_params.sigma_c;
            sigma_s_ind(c) = rf_params.sigma_s;
            diam63_c_ind(c) = rf_params.diam63_c;
            diam63_s_ind(c) = rf_params.diam63_s;
            
            if c <= 3  % Print first few for debugging
                fprintf('    Cell %d: σ_c=%.1f μm, σ_s=%.1f μm, d63_c=%.1f μm, d63_s=%.1f μm\n', ...
                    c, sigma_c_ind(c), sigma_s_ind(c), diam63_c_ind(c), diam63_s_ind(c));
            end
            
        catch ME
            fprintf('    Cell %d: RF estimation failed - %s\n', c, ME.message);
        end
    end
    
    % Summary statistics for individual cells
    valid_idx = individual_results.valid_cells;
    n_valid = sum(valid_idx);
    
    fprintf('  Individual cell analysis complete: %d/%d cells analyzed successfully\n', n_valid, nCells);
    
    if n_valid > 0
        individual_results.sigma_c_mean = mean(sigma_c_ind(valid_idx), 'omitnan');
        individual_results.sigma_c_sem = std(sigma_c_ind(valid_idx), 'omitnan') / sqrt(n_valid);
        individual_results.sigma_s_mean = mean(sigma_s_ind(valid_idx), 'omitnan');
        individual_results.sigma_s_sem = std(sigma_s_ind(valid_idx), 'omitnan') / sqrt(n_valid);
        individual_results.diam63_c_mean = mean(diam63_c_ind(valid_idx), 'omitnan');
        individual_results.diam63_c_sem = std(diam63_c_ind(valid_idx), 'omitnan') / sqrt(n_valid);
        individual_results.diam63_s_mean = mean(diam63_s_ind(valid_idx), 'omitnan');
        individual_results.diam63_s_sem = std(diam63_s_ind(valid_idx), 'omitnan') / sqrt(n_valid);
        
        fprintf('  Center σ: %.1f ± %.1f μm (mean ± SEM)\n', ...
            individual_results.sigma_c_mean, individual_results.sigma_c_sem);
        fprintf('  Surround σ: %.1f ± %.1f μm (mean ± SEM)\n', ...
            individual_results.sigma_s_mean, individual_results.sigma_s_sem);
    end
    
    % ============== GROUP AVERAGE RF ANALYSIS ==============
    fprintf('  Analyzing group average...\n');
    
    % Calculate group mean and SEM across cells
    group_mean = mean(response_matrix, 2, 'omitnan');
    group_sem = std(response_matrix, 0, 2, 'omitnan') ./ sqrt(sum(~isnan(response_matrix), 2));
    
    % Remove sizes with insufficient data
    min_cells_required = max(3, round(0.3 * nCells));  % At least 30% of cells or 3 cells
    valid_sizes = sum(~isnan(response_matrix), 2) >= min_cells_required;
    
    if sum(valid_sizes) < 4
        fprintf('  Group average: Insufficient data points for RF fitting\n');
        group_rf_params = struct();
    else
        % Prepare data for group RF estimation
        rvals_group = sorted_radii(valid_sizes);
        Rmean_group = group_mean(valid_sizes)';
        Rsem_group = group_sem(valid_sizes)';
        
        try
            % Estimate RF parameters for group average
            group_rf_params = estimate_rf_sizes(rvals_group, Rmean_group, Rsem_group, bigR, draws);
            
            fprintf('  Group average: σ_c=%.1f μm (CI: %.1f-%.1f), σ_s=%.1f μm (CI: %.1f-%.1f)\n', ...
                group_rf_params.sigma_c, group_rf_params.ci_c(1), group_rf_params.ci_c(2), ...
                group_rf_params.sigma_s, group_rf_params.ci_s(1), group_rf_params.ci_s(2));
            
        catch ME
            fprintf('  Group average: RF estimation failed - %s\n', ME.message);
            group_rf_params = struct();
        end
    end
    
    % ============== STORE RESULTS ==============
    rf_results.(gname) = struct();
    rf_results.(gname).individual = individual_results;
    rf_results.(gname).group_average = group_rf_params;
    rf_results.(gname).group_mean_response = group_mean;
    rf_results.(gname).group_sem_response = group_sem;
    rf_results.(gname).spot_radii = sorted_radii;
    rf_results.(gname).response_matrix = response_matrix;
    rf_results.(gname).valid_sizes = valid_sizes;
    
    % Collect data for cross-group analysis
    if n_valid > 0
        for idx = find(valid_idx)'
            all_individual_results(end+1, :) = [g, sigma_c_ind(idx), sigma_s_ind(idx), ...
                diam63_c_ind(idx), diam63_s_ind(idx)]; %#ok<AGROW>
            all_group_labels{end+1} = gname; %#ok<AGROW>
        end
    end
end

% ================== VISUALIZATION ==================
fprintf('\n========== CREATING VISUALIZATIONS ==========\n');

% ============== Plot 1: Individual Cell RF Parameters Bar Charts ==============
if ~isempty(all_individual_results)
    figure('Color', 'w', 'Position', [100, 100, 1400, 800]);
    
    % Prepare data for plotting
    group_positions = 1:numel(groups);
    center_means = nan(numel(groups), 1);
    center_sems = nan(numel(groups), 1);
    surround_means = nan(numel(groups), 1);
    surround_sems = nan(numel(groups), 1);
    
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(rf_results, gname) && isfield(rf_results.(gname).individual, 'sigma_c_mean')
            center_means(g) = rf_results.(gname).individual.sigma_c_mean;
            center_sems(g) = rf_results.(gname).individual.sigma_c_sem;
            surround_means(g) = rf_results.(gname).individual.sigma_s_mean;
            surround_sems(g) = rf_results.(gname).individual.sigma_s_sem;
        end
    end
    
    % Center sigma subplot
    subplot(2, 2, 1);
    hold on;
    bar_width = 0.6;
    for g = 1:numel(groups)
        if ~isnan(center_means(g)) && group_colors.isKey(groups{g})
            bar(g, center_means(g), bar_width, 'FaceColor', group_colors(groups{g}), ...
                'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
    errorbar(group_positions, center_means, center_sems, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'CapSize', 8);
    
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    xtickangle(45);
    ylabel('Center σ (μm)');
    title('RF Center Sigma (Individual Cells)');
    grid on;
    
    % Surround sigma subplot  
    subplot(2, 2, 2);
    hold on;
    for g = 1:numel(groups)
        if ~isnan(surround_means(g)) && group_colors.isKey(groups{g})
            bar(g, surround_means(g), bar_width, 'FaceColor', group_colors(groups{g}), ...
                'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
    errorbar(group_positions, surround_means, surround_sems, 'k', 'LineStyle', 'none', ...
        'LineWidth', 1.5, 'CapSize', 8);
    
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    xtickangle(45);
    ylabel('Surround σ (μm)');
    title('RF Surround Sigma (Individual Cells)');
    grid on;
    
    % Center/Surround ratio
    cs_ratio = center_means ./ surround_means;
    subplot(2, 2, 3);
    hold on;
    for g = 1:numel(groups)
        if ~isnan(cs_ratio(g)) && group_colors.isKey(groups{g})
            bar(g, cs_ratio(g), bar_width, 'FaceColor', group_colors(groups{g}), ...
                'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
    
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    xtickangle(45);
    ylabel('Center/Surround σ Ratio');
    title('RF Center/Surround Ratio');
    grid on;
    
    % Sample sizes
    subplot(2, 2, 4);
    n_cells_valid = nan(numel(groups), 1);
    for g = 1:numel(groups)
        gname = groups{g};
        if isfield(rf_results, gname)
            n_cells_valid(g) = sum(rf_results.(gname).individual.valid_cells);
        end
    end
    
    hold on;
    for g = 1:numel(groups)
        if ~isnan(n_cells_valid(g)) && group_colors.isKey(groups{g})
            bar(g, n_cells_valid(g), bar_width, 'FaceColor', group_colors(groups{g}), ...
                'EdgeColor', 'k', 'LineWidth', 1);
        end
    end
    
    xlim([0.5, numel(groups) + 0.5]);
    xticks(group_positions);
    xticklabels(strrep(groups, '_', '\_'));
    xtickangle(45);
    ylabel('Number of Valid Cells');
    title('Sample Sizes');
    grid on;
    
    sgtitle(sprintf('RF Parameters from Individual Cell Analysis (%s)', test_type), ...
        'FontSize', 14, 'FontWeight', 'bold');
    
    % Save figure
    saveas(gcf, fullfile(save_folder, sprintf('RF_IndividualCells_Summary_%s.png', test_type)));
    print(gcf, fullfile(save_folder, sprintf('RF_IndividualCells_Summary_%s', test_type)), '-depsc', '-painters');
end

% ============== Plot 2: Group Average RF Parameters ==============
figure('Color', 'w', 'Position', [200, 150, 1400, 600]);

% Extract group average parameters
group_center_sigma = nan(numel(groups), 1);
group_center_ci_low = nan(numel(groups), 1);
group_center_ci_high = nan(numel(groups), 1);
group_surround_sigma = nan(numel(groups), 1);
group_surround_ci_low = nan(numel(groups), 1);
group_surround_ci_high = nan(numel(groups), 1);

for g = 1:numel(groups)
    gname = groups{g};
    if isfield(rf_results, gname) && ~isempty(fieldnames(rf_results.(gname).group_average))
        group_center_sigma(g) = rf_results.(gname).group_average.sigma_c;
        group_center_ci_low(g) = rf_results.(gname).group_average.ci_c(1);
        group_center_ci_high(g) = rf_results.(gname).group_average.ci_c(2);
        group_surround_sigma(g) = rf_results.(gname).group_average.sigma_s;
        group_surround_ci_low(g) = rf_results.(gname).group_average.ci_s(1);
        group_surround_ci_high(g) = rf_results.(gname).group_average.ci_s(2);
    end
end

% Center sigma
subplot(1, 3, 1);
hold on;
for g = 1:numel(groups)
    if ~isnan(group_center_sigma(g)) && group_colors.isKey(groups{g})
        bar(g, group_center_sigma(g), 0.6, 'FaceColor', group_colors(groups{g}), ...
            'EdgeColor', 'k', 'LineWidth', 1);
        
        % Add confidence intervals
        errorbar(g, group_center_sigma(g), ...
            group_center_sigma(g) - group_center_ci_low(g), ...
            group_center_ci_high(g) - group_center_sigma(g), ...
            'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
    end
end

xlim([0.5, numel(groups) + 0.5]);
xticks(1:numel(groups));
xticklabels(strrep(groups, '_', '\_'));
xtickangle(45);
ylabel('Center σ (μm)');
title('RF Center Sigma (Group Average)');
grid on;

% Surround sigma
subplot(1, 3, 2);
hold on;
for g = 1:numel(groups)
    if ~isnan(group_surround_sigma(g)) && group_colors.isKey(groups{g})
        bar(g, group_surround_sigma(g), 0.6, 'FaceColor', group_colors(groups{g}), ...
            'EdgeColor', 'k', 'LineWidth', 1);
        
        % Add confidence intervals
        errorbar(g, group_surround_sigma(g), ...
            group_surround_sigma(g) - group_surround_ci_low(g), ...
            group_surround_ci_high(g) - group_surround_sigma(g), ...
            'k', 'LineStyle', 'none', 'LineWidth', 1.5, 'CapSize', 8);
    end
end

xlim([0.5, numel(groups) + 0.5]);
xticks(1:numel(groups));
xticklabels(strrep(groups, '_', '\_'));
xtickangle(45);
ylabel('Surround σ (μm)');
title('RF Surround Sigma (Group Average)');
grid on;

% Center/Surround ratio
subplot(1, 3, 3);
group_cs_ratio = group_center_sigma ./ group_surround_sigma;
hold on;
for g = 1:numel(groups)
    if ~isnan(group_cs_ratio(g)) && group_colors.isKey(groups{g})
        bar(g, group_cs_ratio(g), 0.6, 'FaceColor', group_colors(groups{g}), ...
            'EdgeColor', 'k', 'LineWidth', 1);
    end
end

xlim([0.5, numel(groups) + 0.5]);
xticks(1:numel(groups));
xticklabels(strrep(groups, '_', '\_'));
xtickangle(45);
ylabel('Center/Surround σ Ratio');
title('RF Center/Surround Ratio (Group Average)');
grid on;

sgtitle(sprintf('RF Parameters from Group Average Analysis (%s)', test_type), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, fullfile(save_folder, sprintf('RF_GroupAverage_Summary_%s.png', test_type)));
print(gcf, fullfile(save_folder, sprintf('RF_GroupAverage_Summary_%s', test_type)), '-depsc', '-painters');

% ============== Plot 3: Response Curves with RF Fits ==============
figure('Color', 'w', 'Position', [300, 200, 1600, 800]);

for g = 1:numel(groups)
    gname = groups{g};
    if ~isfield(rf_results, gname)
        continue;
    end
    
    subplot(2, numel(groups), g);
    
    % Plot individual cell responses (light gray)
    response_matrix = rf_results.(gname).response_matrix;
    hold on;
    for c = 1:size(response_matrix, 2)
        plot(sorted_radii, response_matrix(:, c), '-', 'Color', [0.8, 0.8, 0.8], 'LineWidth', 0.5);
    end
    
    % Plot group mean ± SEM
    group_mean = rf_results.(gname).group_mean_response;
    group_sem = rf_results.(gname).group_sem_response;
    
    if group_colors.isKey(gname)
        color = group_colors(gname);
    else
        color = [0, 0, 0];
    end
    
    % Shaded error region
    valid_idx = ~isnan(group_mean) & ~isnan(group_sem);
    if sum(valid_idx) > 0
        fill([sorted_radii(valid_idx), fliplr(sorted_radii(valid_idx))], ...
             [group_mean(valid_idx)' + group_sem(valid_idx)', ...
              fliplr(group_mean(valid_idx)' - group_sem(valid_idx)')], ...
             color, 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    end
    
    % Group mean line
    plot(sorted_radii, group_mean, '-', 'Color', color, 'LineWidth', 2);
    
    xlabel('Spot Radius (μm)');
    ylabel('Response (spikes/s)');
    title(strrep(gname, '_', '\_'));
    grid on;
    xlim([0, max(sorted_radii)]);
    
    % Add RF parameter text if available
    if ~isempty(fieldnames(rf_results.(gname).group_average))
        rf_params = rf_results.(gname).group_average;
        text(0.05, 0.95, sprintf('σ_c = %.1f μm\nσ_s = %.1f μm', ...
            rf_params.sigma_c, rf_params.sigma_s), ...
            'Units', 'normalized', 'VerticalAlignment', 'top', ...
            'BackgroundColor', [1 1 1 0.8], 'FontSize', 9);
    end
    
    % Lower subplot: Individual cell distribution
    subplot(2, numel(groups), g + numel(groups));
    
    individual_data = rf_results.(gname).individual;
    if sum(individual_data.valid_cells) > 0
        valid_cells = find(individual_data.valid_cells);
        sigma_c_values = nan(length(valid_cells), 1);
        sigma_s_values = nan(length(valid_cells), 1);
        
        for i = 1:length(valid_cells)
            c = valid_cells(i);
            if ~isempty(individual_data.rf_params{c})
                sigma_c_values(i) = individual_data.rf_params{c}.sigma_c;
                sigma_s_values(i) = individual_data.rf_params{c}.sigma_s;
            end
        end
        
        scatter(sigma_c_values, sigma_s_values, 50, color, 'filled', 'MarkerEdgeColor', 'k');
        xlabel('Center σ (μm)');
        ylabel('Surround σ (μm)');
        title(sprintf('Individual Cells (n=%d)', sum(~isnan(sigma_c_values))));
        grid on;
        
        % Add unity line
        max_val = max([sigma_c_values; sigma_s_values]);
        plot([0, max_val], [0, max_val], 'k--', 'LineWidth', 1);
    else
        text(0.5, 0.5, 'No valid cells', 'Units', 'normalized', ...
            'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle');
    end
end

sgtitle(sprintf('Response Curves and RF Parameter Distributions (%s)', test_type), ...
    'FontSize', 14, 'FontWeight', 'bold');

% Save figure
saveas(gcf, fullfile(save_folder, sprintf('RF_ResponseCurves_%s.png', test_type)));

% ================== SAVE RESULTS ==================
fprintf('\nSaving results...\n');

% Save comprehensive results structure
save(fullfile(save_folder, sprintf('RF_Analysis_Results_%s.mat', test_type)), ...
    'rf_results', 'groups', 'test_type', 'draws', 'bigR', 'sorted_radii', 'sorted_diameters');

% ================== SUMMARY STATISTICS ==================
fprintf('\n========== RF ANALYSIS SUMMARY (%s) ==========\n', test_type);

for g = 1:numel(groups)
    gname = groups{g};
    if ~isfield(rf_results, gname)
        continue;
    end
    
    fprintf('\n%s:\n', gname);
    
    % Individual cell results
    individual_data = rf_results.(gname).individual;
    n_valid = sum(individual_data.valid_cells);
    n_total = individual_data.n_cells;
    
    fprintf('  Individual cells: %d/%d analyzed successfully (%.1f%%)\n', ...
        n_valid, n_total, 100*n_valid/n_total);
    
    if n_valid > 0
        fprintf('    Center σ: %.1f ± %.1f μm (mean ± SEM)\n', ...
            individual_data.sigma_c_mean, individual_data.sigma_c_sem);
        fprintf('    Surround σ: %.1f ± %.1f μm (mean ± SEM)\n', ...
            individual_data.sigma_s_mean, individual_data.sigma_s_sem);
        fprintf('    Center diameter (63%% mass): %.1f ± %.1f μm\n', ...
            individual_data.diam63_c_mean, individual_data.diam63_c_sem);
        fprintf('    Surround diameter (63%% mass): %.1f ± %.1f μm\n', ...
            individual_data.diam63_s_mean, individual_data.diam63_s_sem);
    end
    
    % Group average results
    if ~isempty(fieldnames(rf_results.(gname).group_average))
        group_params = rf_results.(gname).group_average;
        fprintf('  Group average:\n');
        fprintf('    Center σ: %.1f μm (68%% CI: %.1f-%.1f)\n', ...
            group_params.sigma_c, group_params.ci_c(1), group_params.ci_c(2));
        fprintf('    Surround σ: %.1f μm (68%% CI: %.1f-%.1f)\n', ...
            group_params.sigma_s, group_params.ci_s(1), group_params.ci_s(2));
        fprintf('    Center diameter (63%% mass): %.1f μm\n', group_params.diam63_c);
        fprintf('    Surround diameter (63%% mass): %.1f μm\n', group_params.diam63_s);
    else
        fprintf('  Group average: RF fitting failed\n');
    end
end

fprintf('\nAnalysis complete! Results saved to: %s\n', save_folder);
fprintf('Generated files:\n');
fprintf('  - RF_Analysis_Results_%s.mat (complete results)\n', test_type);
fprintf('  - RF_IndividualCells_Summary_%s.png/eps (individual cell summary)\n', test_type);
fprintf('  - RF_GroupAverage_Summary_%s.png/eps (group average summary)\n', test_type);
fprintf('  - RF_ResponseCurves_%s.png (response curves and distributions)\n', test_type);
