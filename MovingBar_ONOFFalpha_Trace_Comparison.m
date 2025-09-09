close all; clear; clc;
%% Load processed data
set_name = 'latest';  % before081425
is_show_fitted = 1;
switch set_name
     case 'before081425'
        data_sets = {'e100724', 'f100724', 'a101224', 'b101224', 'c101224',   'd101224', 'e101224',...
                    'b101424', 'c101424', 'd101424', 'e101424', 'a101624',   'b101624', 'd101624', 'e101624',...
                    'b101924', 'c101924', 'd101924', 'e101924'};
        cell_type = {'OFF',      'OFF',    'OFF',      'ON',       'OFF',     'ON',      'OFF',...
                    'OFF',      'OFF',    'ON',       'ON',       'ON',      'ON',      'ON',       'ON',...
                    'ON',       'OFF',    'OFF',      'OFF'};
        location =  {'Temporal', 'Temporal','Nasal',   'Nasal',    'Nasal',   'Nasal',   'Nasal',...
                    'Temporal', 'Temporal','Temporal','Temporal', 'Nasal',   'Nasal',   'Nasal',    'Nasal',...
                    'Nasal',    'Nasal',   'Nasal',   'Nasal'};

    case 'latest'
        data_sets = {'e100724', 'f100724', 'a101224', 'b101224', 'c101224',   'd101224', 'e101224',...
                    'b101424', 'c101424', 'd101424', 'e101424', 'a101624',   'b101624', 'd101624', 'e101624',...
                    'b101924', 'c101924', 'd101924', 'e101924', 'b103124',   'e103124', 'a110424',...
                    'c110424',                       'f110424', 'g110424',   'a110924', 'b110924', 'c110924',...
                    'a111224'};
        cell_type = {'OFF',      'OFF',    'OFF',      'ON',       'OFF',     'ON',      'OFF',...
                    'OFF',      'OFF',    'ON',       'ON',       'ON',      'ON',      'ON',       'ON',...
                    'ON',       'OFF',    'OFF',      'OFF',      'ON',      'OFF',     'ON',...
                    'ON',                             'ON',       'OFF',     'ON',      'OFF',      'OFF',...
                    'ON'};
        location =  {'Temporal', 'Temporal','Nasal',   'Nasal',    'Nasal',   'Nasal',   'Nasal',...
                    'Temporal', 'Temporal','Temporal','Temporal', 'Nasal',   'Nasal',   'Nasal',    'Nasal',...
                    'Nasal',    'Nasal',   'Nasal',   'Nasal',    'Temporal','Temporal','Temporal', ...
                    'Temporal',                      'Temporal',  'Temporal','Temporal','Temporal', 'Temporal',...
                    'Temporal'};
end

save_fig_folder = './Figures/MovingBarSummary/';
if ~exist(save_fig_folder, 'dir')
    mkdir(save_fig_folder);
end
clear Data 
num_set = length(data_sets);
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar';
Cdat = nan(num_set, 10);
Cbas = nan(num_set, 1);
Csw = nan(num_set, 4);
% Add LNK_params array to collect parameters
LNK_params_array = nan(num_set, 10); % 10 columns for the specified parameters
implement_case_id = 6;
for i = 1:num_set
    file_name = sprintf('%s_moving_bar_processed.mat', data_sets{i});
    Data{i} = load(fullfile(folder_name, file_name));
    if is_show_fitted
        file_name = sprintf('%s_moving_bar_fitted.mat', data_sets{i});
        loaded_data = load(sprintf('./Results/MovingBar/%s', file_name), 'PredictionResults',...
         'BaselineCorr', 'LNK_params_s', 'LNK_params_w', 'LN_params_s', 'LNK_params_d');
        Cdat(i, :) = loaded_data.PredictionResults;
        Cbas(i, :) = loaded_data.BaselineCorr;
        Csw(i, :) = [loaded_data.LNK_params_s.w_xs loaded_data.LNK_params_w.w_xs,...
                     loaded_data.LN_params_s.gamma loaded_data.LNK_params_d.w_xs];

        % Collect LNK_params_w parameters in specified order
        % Order: 'tau', 'alpha_d', 'theta', 'sigma0', 'alpha', 'beta', 'b_out', 'g_out', 'w_xs', 'dt'
        LNK_params_array(i, 1) = getfield_safe(loaded_data.LNK_params_w, 'tau');
        LNK_params_array(i, 2) = getfield_safe(loaded_data.LNK_params_w, 'alpha_d');
        LNK_params_array(i, 3) = getfield_safe(loaded_data.LNK_params_w, 'theta');
        LNK_params_array(i, 4) = getfield_safe(loaded_data.LNK_params_w, 'sigma0');
        LNK_params_array(i, 5) = getfield_safe(loaded_data.LNK_params_w, 'alpha');
        LNK_params_array(i, 6) = getfield_safe(loaded_data.LNK_params_w, 'beta');
        LNK_params_array(i, 7) = getfield_safe(loaded_data.LNK_params_w, 'b_out');
        LNK_params_array(i, 8) = getfield_safe(loaded_data.LNK_params_w, 'g_out');
        LNK_params_array(i, 9) = getfield_safe(loaded_data.LNK_params_w, 'w_xs');
        LNK_params_array(i, 10) = 0.01; % dt default value
        % Helper function to safely get field or assign NaN
        
    end
end

%% Load fitted parameters from White Noise for LN model
process_version = 'GaussianFitting_processed_082025_1.mat';
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
processedFile = fullfile(folder_name, process_version);

if exist(processedFile, 'file')
    WN = load(processedFile, 'gauss_est', 'Gauss_TF_est', 'data_sets', 'cell_type', 'location', ...
    'gauss_est_ON_nasal', 'gauss_est_ON_temporal', 'Gauss_TF_est_OFF_nasal', 'Gauss_TF_est_OFF_temporal');
    fprintf('Loaded fitted parameters for LN model\n');
else
    error('Fitted parameters not found. Run WhiteNoise_ONOFFalpha_Comparison.m first.');
end

% Check session for consistency between loaded variables and WN
for i = 1:num_set
    ds_name = data_sets{i};
    wn_idx = find(strcmp(ds_name, WN.data_sets));
    if isempty(wn_idx)
        error('WARNING: %s not found in WN.data_sets\n', ds_name);
        continue;
    end
    % If WN contains cell_type and location fields, compare them
    if isfield(WN, 'cell_type')
        if ~strcmp(cell_type{i}, WN.cell_type{wn_idx})
            fprintf('Mismatch in cell_type for %s: %s (main) vs %s (WN)\n', ds_name, cell_type{i}, WN.cell_type{wn_idx});
        end
    else
        error('WN data does not contain cell_type field');
    end
    if isfield(WN, 'location')
        if ~strcmp(location{i}, WN.location{wn_idx})
            fprintf('Mismatch in location for %s: %s (main) vs %s (WN)\n', ds_name, location{i}, WN.location{wn_idx});
        end
    else
        error('WN data does not contain location field');
    end
end
fprintf('Check session completed.\n');

%% Split LNK parameters into groups based on WN fields
if is_show_fitted && exist('WN', 'var')
    % Create mapping from data_sets to WN indices
    LNK_params_groups = struct();
    
    % Initialize arrays for each group
    LNK_params_ON_temporal = [];
    LNK_params_ON_nasal = [];
    LNK_params_OFF_temporal = [];
    LNK_params_OFF_nasal = [];
    
    % Iterate through WN.data_sets to ensure consistent ordering and size
    for wn_idx = 1:length(WN.data_sets)
        wn_ds_name = WN.data_sets{wn_idx};
        wn_cell_type = WN.cell_type{wn_idx};
        wn_location = WN.location{wn_idx};
        
        % Find corresponding index in main data_sets
        main_idx = find(strcmp(wn_ds_name, data_sets));
        
        if isempty(main_idx)
            % Dataset exists in WN but not in main script - use NaN row
            fprintf('INFO: %s found in WN but not in main data_sets, using NaN values\n', wn_ds_name);
            param_row = nan(1, 10); % 10 parameters with NaN values
        else
            % Dataset exists in both - use actual parameters
            param_row = LNK_params_array(main_idx, :);
        end
        
        % Assign to appropriate group based on WN cell type and location
        if strcmpi(wn_cell_type, 'ON') && strcmpi(wn_location, 'Temporal')
            LNK_params_ON_temporal = [LNK_params_ON_temporal; param_row];
        elseif strcmpi(wn_cell_type, 'ON') && strcmpi(wn_location, 'Nasal')
            LNK_params_ON_nasal = [LNK_params_ON_nasal; param_row];
        elseif strcmpi(wn_cell_type, 'OFF') && strcmpi(wn_location, 'Temporal')
            LNK_params_OFF_temporal = [LNK_params_OFF_temporal; param_row];
        elseif strcmpi(wn_cell_type, 'OFF') && strcmpi(wn_location, 'Nasal')
            LNK_params_OFF_nasal = [LNK_params_OFF_nasal; param_row];
        else
            fprintf('WARNING: Unknown cell type/location combination: %s/%s for %s\n', wn_cell_type, wn_location, wn_ds_name);
        end
    end
    
    % Store in structure with parameter names
    param_names = {'tau', 'alpha_d', 'theta', 'sigma0', 'alpha', 'beta', 'b_out', 'g_out', 'w_xs', 'dt'};
    LNK_params_groups.param_names = param_names;
    LNK_params_groups.ON_temporal = LNK_params_ON_temporal;
    LNK_params_groups.ON_nasal = LNK_params_ON_nasal;
    LNK_params_groups.OFF_temporal = LNK_params_OFF_temporal;
    LNK_params_groups.OFF_nasal = LNK_params_OFF_nasal;
    
    % Sanity check with WN group sizes
    fprintf('\n=== LNK PARAMETERS SPLIT SANITY CHECK ===\n');
    
    % Check if WN group variables exist and compare sizes
    if isfield(WN, 'gauss_est_ON_nasal')
        expected_ON_nasal = size(WN.gauss_est_ON_nasal, 1);
        actual_ON_nasal = size(LNK_params_ON_nasal, 1);
        fprintf('ON-nasal: Expected %d, Got %d %s\n', expected_ON_nasal, actual_ON_nasal, ...
                iif(expected_ON_nasal == actual_ON_nasal, '✓', '✗'));
    end
    
    if isfield(WN, 'gauss_est_ON_temporal')
        expected_ON_temporal = size(WN.gauss_est_ON_temporal, 1);
        actual_ON_temporal = size(LNK_params_ON_temporal, 1);
        fprintf('ON-temporal: Expected %d, Got %d %s\n', expected_ON_temporal, actual_ON_temporal, ...
                iif(expected_ON_temporal == actual_ON_temporal, '✓', '✗'));
    end
    
    if isfield(WN, 'Gauss_TF_est_OFF_nasal')
        expected_OFF_nasal = size(WN.Gauss_TF_est_OFF_nasal, 1);
        actual_OFF_nasal = size(LNK_params_OFF_nasal, 1);
        fprintf('OFF-nasal: Expected %d, Got %d %s\n', expected_OFF_nasal, actual_OFF_nasal, ...
                iif(expected_OFF_nasal == actual_OFF_nasal, '✓', '✗'));
    end
    
    if isfield(WN, 'Gauss_TF_est_OFF_temporal')
        expected_OFF_temporal = size(WN.Gauss_TF_est_OFF_temporal, 1);
        actual_OFF_temporal = size(LNK_params_OFF_temporal, 1);
        fprintf('OFF-temporal: Expected %d, Got %d %s\n', expected_OFF_temporal, actual_OFF_temporal, ...
                iif(expected_OFF_temporal == actual_OFF_temporal, '✓', '✗'));
    end
    
    fprintf('\nLNK parameters successfully split into groups.\n');
    fprintf('Use LNK_params_groups.ON_temporal, LNK_params_groups.ON_nasal, etc. to access the data.\n');
    fprintf('Column order: %s\n', strjoin(param_names, ', '));
    
else
    fprintf('LNK parameter splitting skipped: is_show_fitted=%d, WN exists=%d\n', is_show_fitted, exist('WN', 'var'));
end

cell_type_numeric = cellfun(@(x) strcmp(x, 'ON'), cell_type);
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), location);

keyboard;
%% Plot Cdat (4th and 5th columns) and Cbas
if is_show_fitted
    % Define groups
    groups = {'ON-Temporal', 'ON-Nasal', 'OFF-Temporal', 'OFF-Nasal'};
    group_indices = {
        cell_type_numeric == 1 & location_type_numeric == 1; ... % ON-Temporal
        cell_type_numeric == 1 & location_type_numeric == 0; ... % ON-Nasal
        cell_type_numeric == 0 & location_type_numeric == 1; ... % OFF-Temporal
        cell_type_numeric == 0 & location_type_numeric == 0  ... % OFF-Nasal
    };
    
    % Figure for Cdat and Cbas
    figure('Name', 'Model Performance Comparison');
    for g = 1:4
        subplot(2, 2, g);
        hold on;
        
        % Get data for this group
        group_mask = group_indices{g};
        if sum(group_mask) == 0
            title([groups{g} ' (n=0)']);
            continue;
        end
        
        group_cbas = Cbas(group_mask);
        group_cdat_ln = Cdat(group_mask, 5);    % 5th column: LN
        group_cdat_lnk = Cdat(group_mask, 9);   % 9th column: LNK

        % Sort by Cbas in descending order
        [sorted_cbas, sort_idx] = sort(group_cbas, 'descend');
        sorted_cdat_ln = group_cdat_ln(sort_idx);
        sorted_cdat_lnk = group_cdat_lnk(sort_idx);
        
        x_vals = 1:length(sorted_cbas);
        
        % Plot lines
        plot(x_vals, sorted_cbas, 'k-', 'LineWidth', 2, 'DisplayName', 'Repeat reliability');
        plot(x_vals, sorted_cdat_ln, 'b-', 'LineWidth', 2, 'DisplayName', 'LN');
        plot(x_vals, sorted_cdat_lnk, 'r-', 'LineWidth', 2, 'DisplayName', 'LNK');
        
        xlabel('Cell (sorted by repeat reliability)');
        ylabel('Correlation coefficient');
        title([groups{g} ' (n=' num2str(sum(group_mask)) ')']);
        xticks(1:2:length(sorted_cbas));
        xticklabels(arrayfun(@num2str, 1:2:length(sorted_cbas), 'UniformOutput', false));
        yticks(0:0.25:1);
        yticklabels(arrayfun(@num2str, 0:0.25:1, 'UniformOutput', false));
        ylim([0 1]);
        legend('Location', 'best');
        grid off;
    end
    sgtitle('Model Performance: LN vs LNK vs Repeat Reliability');

    saveas(gcf, fullfile(save_fig_folder, 'Model-Performance-LN-LNK.png'));
    print(gcf, fullfile(save_fig_folder, 'Model-Performance-LN-LNK'), '-depsc', '-painters');

    % Figure for Csw (bar chart) - 2 subplots for ON/OFF with 2 bars each (Temporal vs Nasal averages)
    figure('Name', 'Subunit Weights (Csw)');
    
    % Define cell types and locations
    cell_types = {'ON', 'OFF'};
    locations = {'Temporal', 'Nasal'};
    colors = [0.4 0.6 0.8; 0.8 0.4 0.6]; % Colors for Temporal and Nasal
    
    for ct = 1:2 % Loop through ON and OFF
        subplot(1, 2, ct);
        
        % Get data for each location within this cell type
        temporal_mask = (cell_type_numeric == (ct==1)) & (location_type_numeric == 1);
        nasal_mask = (cell_type_numeric == (ct==1)) & (location_type_numeric == 0);
        
        temporal_csw = Csw(temporal_mask, 3);
        nasal_csw = Csw(nasal_mask, 3);
        
        % Remove outliers using interquartile range method
        if length(temporal_csw) > 0
            Q1_t = prctile(temporal_csw, 25);
            Q3_t = prctile(temporal_csw, 75);
            IQR_t = Q3_t - Q1_t;
            temporal_csw_clean = temporal_csw(temporal_csw >= Q1_t - 1.5*IQR_t & temporal_csw <= Q3_t + 1.5*IQR_t);
        else
            temporal_csw_clean = [];
        end
        
        if length(nasal_csw) > 0
            Q1_n = prctile(nasal_csw, 25);
            Q3_n = prctile(nasal_csw, 75);
            IQR_n = Q3_n - Q1_n;
            nasal_csw_clean = nasal_csw(nasal_csw >= Q1_n - 1.5*IQR_n & nasal_csw <= Q3_n + 1.5*IQR_n);
        else
            nasal_csw_clean = [];
        end
        
        % Calculate means
        mean_temporal = mean(temporal_csw_clean);
        mean_nasal = mean(nasal_csw_clean);
        
        % Calculate standard errors
        se_temporal = std(temporal_csw_clean) / sqrt(length(temporal_csw_clean));
        se_nasal = std(nasal_csw_clean) / sqrt(length(nasal_csw_clean));
        
        if isnan(mean_temporal), mean_temporal = 0; se_temporal = 0; end
        if isnan(mean_nasal), mean_nasal = 0; se_nasal = 0; end
        
        % Create bar chart with only 2 bars
        bar_data = [mean_temporal, mean_nasal];
        error_data = [se_temporal, se_nasal];
        
        bars = bar(bar_data);
        bars.FaceColor = 'flat';
        bars.CData(1,:) = colors(1, :); % Temporal
        bars.CData(2,:) = colors(2, :); % Nasal
        
        hold on;
        
        % Add error bars
        errorbar(1:2, bar_data, error_data, 'k', 'LineStyle', 'none', 'LineWidth', 1.5);
        
        % Add individual data points as scatter
        x_temporal = ones(size(temporal_csw_clean)) + 0.1*randn(size(temporal_csw_clean))*0.1; % Add small jitter
        x_nasal = 2*ones(size(nasal_csw_clean)) + 0.1*randn(size(nasal_csw_clean))*0.1;
        
        scatter(x_temporal, temporal_csw_clean, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
        scatter(x_nasal, nasal_csw_clean, 30, 'k', 'filled', 'MarkerFaceAlpha', 0.6);
        
        % Formatting
        set(gca, 'XTick', 1:2, 'XTickLabel', locations);
        ylabel('Subunit Weight');
        title([cell_types{ct} ' (Temporal: n=' num2str(length(temporal_csw_clean)) ...
               '/' num2str(sum(temporal_mask)) ', Nasal: n=' num2str(length(nasal_csw_clean)) ...
               '/' num2str(sum(nasal_mask)) ')']);
        grid on;
        
        % Set y-axis limits to show data nicely
        all_data = [temporal_csw_clean; nasal_csw_clean];
        % if ~isempty(all_data)
        %     ylim([0, max(all_data)*1.2]);
        % end
    end
    sgtitle('Average Subunit Weights by Cell Type and Location (Outliers Removed)');
    
    %% Statistical Testing for Subunit Weights
    fprintf('\n=== STATISTICAL ANALYSIS OF SUBUNIT WEIGHTS ===\n');
    
    % Prepare data for statistical testing
    stat_results = struct();
    
    for ct = 1:2 % Loop through ON and OFF
        cell_type_name = cell_types{ct};
        fprintf('\n--- %s CELLS ---\n', upper(cell_type_name));
        
        % Get data for each location within this cell type
        temporal_mask = (cell_type_numeric == (ct==1)) & (location_type_numeric == 1);
        nasal_mask = (cell_type_numeric == (ct==1)) & (location_type_numeric == 0);
        
        temporal_csw = Csw(temporal_mask, 2);
        nasal_csw = Csw(nasal_mask, 2);
        
        % Remove outliers using same method as in plotting
        if length(temporal_csw) > 0
            Q1_t = prctile(temporal_csw, 25);
            Q3_t = prctile(temporal_csw, 75);
            IQR_t = Q3_t - Q1_t;
            temporal_csw_clean = temporal_csw(temporal_csw >= Q1_t - 1.5*IQR_t & temporal_csw <= Q3_t + 1.5*IQR_t);
        else
            temporal_csw_clean = [];
        end
        
        if length(nasal_csw) > 0
            Q1_n = prctile(nasal_csw, 25);
            Q3_n = prctile(nasal_csw, 75);
            IQR_n = Q3_n - Q1_n;
            nasal_csw_clean = nasal_csw(nasal_csw >= Q1_n - 1.5*IQR_n & nasal_csw <= Q3_n + 1.5*IQR_n);
        else
            nasal_csw_clean = [];
        end
        
        % Store cleaned data
        stat_results.(sprintf('%s_temporal', lower(cell_type_name))) = temporal_csw_clean;
        stat_results.(sprintf('%s_nasal', lower(cell_type_name))) = nasal_csw_clean;
        
        % Descriptive statistics
        fprintf('Temporal: n=%d, mean=%.4f, std=%.4f, median=%.4f\n', ...
            length(temporal_csw_clean), mean(temporal_csw_clean), std(temporal_csw_clean), median(temporal_csw_clean));
        fprintf('Nasal: n=%d, mean=%.4f, std=%.4f, median=%.4f\n', ...
            length(nasal_csw_clean), mean(nasal_csw_clean), std(nasal_csw_clean), median(nasal_csw_clean));
        
        % Statistical tests (only if both groups have data)
        if length(temporal_csw_clean) >= 3 && length(nasal_csw_clean) >= 3
            
            % 1. Normality tests (Shapiro-Wilk if available, otherwise Lilliefors)
            try
                % Test normality for each group
                [h_norm_t, p_norm_t] = lillietest(temporal_csw_clean);
                [h_norm_n, p_norm_n] = lillietest(nasal_csw_clean);
                
                fprintf('Normality tests (Lilliefors):\n');
                fprintf('  Temporal: p=%.4f (%s)\n', p_norm_t, iif(h_norm_t, 'NOT normal', 'normal'));
                fprintf('  Nasal: p=%.4f (%s)\n', p_norm_n, iif(h_norm_n, 'NOT normal', 'normal'));
                
                % 2. Equal variance test (F-test)
                [h_var, p_var] = vartest2(temporal_csw_clean, nasal_csw_clean);
                fprintf('Equal variance test (F-test): p=%.4f (%s)\n', p_var, iif(h_var, 'unequal variances', 'equal variances'));
                
                % 3. Choose appropriate test based on normality and variance
                if ~h_norm_t && ~h_norm_n && ~h_var
                    % Both normal and equal variances -> two-sample t-test
                    [h_ttest, p_ttest] = ttest2(temporal_csw_clean, nasal_csw_clean);
                    fprintf('Two-sample t-test: p=%.4f (%s)\n', p_ttest, iif(h_ttest, 'SIGNIFICANT', 'not significant'));
                    stat_results.(sprintf('%s_test_used', lower(cell_type_name))) = 'two-sample t-test';
                    stat_results.(sprintf('%s_p_value', lower(cell_type_name))) = p_ttest;
                    
                elseif ~h_norm_t && ~h_norm_n && h_var
                    % Both normal but unequal variances -> Welch's t-test
                    [h_welch, p_welch] = ttest2(temporal_csw_clean, nasal_csw_clean, 'Vartype', 'unequal');
                    fprintf('Welch''s t-test (unequal variances): p=%.4f (%s)\n', p_welch, iif(h_welch, 'SIGNIFICANT', 'not significant'));
                    stat_results.(sprintf('%s_test_used', lower(cell_type_name))) = 'Welch t-test';
                    stat_results.(sprintf('%s_p_value', lower(cell_type_name))) = p_welch;
                    
                else
                    % Non-normal data -> Mann-Whitney U test (Wilcoxon rank-sum)
                    [p_ranksum, h_ranksum] = ranksum(temporal_csw_clean, nasal_csw_clean);
                    fprintf('Mann-Whitney U test (Wilcoxon rank-sum): p=%.4f (%s)\n', p_ranksum, iif(h_ranksum, 'SIGNIFICANT', 'not significant'));
                    stat_results.(sprintf('%s_test_used', lower(cell_type_name))) = 'Mann-Whitney U';
                    stat_results.(sprintf('%s_p_value', lower(cell_type_name))) = p_ranksum;
                end
                
                % 4. Effect size (Cohen's d)
                pooled_std = sqrt(((length(temporal_csw_clean)-1)*var(temporal_csw_clean) + ...
                                 (length(nasal_csw_clean)-1)*var(nasal_csw_clean)) / ...
                                (length(temporal_csw_clean) + length(nasal_csw_clean) - 2));
                cohens_d = (mean(temporal_csw_clean) - mean(nasal_csw_clean)) / pooled_std;
                fprintf('Effect size (Cohen''s d): %.4f (%s)\n', cohens_d, ...
                    iif(abs(cohens_d) < 0.2, 'negligible', ...
                    iif(abs(cohens_d) < 0.5, 'small', ...
                    iif(abs(cohens_d) < 0.8, 'medium', 'large'))));
                
                stat_results.(sprintf('%s_cohens_d', lower(cell_type_name))) = cohens_d;
                
            catch ME
                fprintf('Error in statistical testing: %s\n', ME.message);
            end
            
        else
            fprintf('Insufficient data for statistical testing (need n>=3 for both groups)\n');
        end
    end
    
    % Save statistical results
    save('subunit_weight_statistics.mat', 'stat_results');
    fprintf('\nStatistical results saved to subunit_weight_statistics.mat\n');
end

% Helper function for conditional text
function result = iif(condition, true_text, false_text)
    if condition
        result = true_text;
    else
        result = false_text;
    end
end

%
keyboard;
%%
%%
if is_show_fitted
    is_blurry = 0;
    Colors = lines(4);
    figure; 
    subplot(1, 2, 1);hold on
    cavg = mean(Cbas(cell_type_numeric==1, :), 2);
    plot(cavg, 'k');
    cavg = mean(Cdat(cell_type_numeric==1, 1, 1), 2);
    plot(cavg, 'Color',Colors(1, :));
    if is_blurry
        cavg = mean(Cdat(cell_type_numeric==1, 2, 1), 2);
        plot(cavg, 'Color',Colors(2, :));
        legend({'Repeat reliability', 'Standard', 'Blurry'})
    else
        legend({'Repeat reliability', 'Standard'})
    end
    xlabel('# cell');
    ylabel('Corr. coeff.');
    ylim([0 1])

    title('ON');

    subplot(1, 2, 2);hold on
    cavg = mean(Cbas(cell_type_numeric==0, :), 2);
    plot(cavg, 'k');
    cavg = mean(Cdat(cell_type_numeric==0, 1, 1), 2);
    plot(cavg, 'Color',Colors(1, :));
    if is_blurry
        cavg = mean(Cdat(cell_type_numeric==0, 2, 1), 2);
        plot(cavg, 'Color',Colors(2, :));
        legend({'Repeat reliability', 'Standard', 'Blurry'})
    else
        legend({'Repeat reliability', 'Standard'})
    end
    xlabel('# cell');
    ylabel('Corr. coeff.');
    ylim([0 1])

    title('OFF');
    %%
    keyboard;
end
%% ON & OFF comparison
Fz = 100;
disp_direction = 180;
disp_contrast = 0.33;
disp_bar_witdth = [50, 100, 200, 400, 800];
disp_speeds = [500, 1000, 2000, 4000, 8000];
max_t = 459;
ct = (0:max_t-1)/Fz;
Colors = parula(4);
cell_type_numeric = cellfun(@(x) strcmp(x, 'ON'), cell_type);
for i = 1:length(disp_contrast)
    Trace = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
    figure;
    for j = 1:length(disp_bar_witdth)
        for q = 1:length(disp_speeds)
            subplot(length(disp_bar_witdth), length(disp_speeds), (j-1)*length(disp_speeds)+q); hold on

            for k = 1:num_set
                dir_id = find(Data{k}.dim1_moving_direction == disp_direction);
                ctr_id = find(Data{k}.dim2_contrast == disp_contrast);
                bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j));
                sp_id = find(Data{k}.dim4_speeds == disp_speeds(q));
                csig = squeeze(mean(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
                Trace(j, q, k, :) = csig(1:max_t);
                switch lower(cell_type{k})
                    case 'on'
                        plot(ct, csig(1:max_t), 'Color', [247, 224, 12]/255);
                    case 'off'
                        plot(ct, csig(1:max_t), 'Color', 0.4*ones(1, 3));
                end
            end
            plot(ct, squeeze(mean(Trace(j, q, cell_type_numeric==1, :), 3)), 'Color', [245 182 66]/255, 'LineWidth', 2)
            plot(ct, squeeze(mean(Trace(j, q, cell_type_numeric==0, :), 3)), 'Color', 0*ones(1, 3), 'LineWidth', 2)
            ylim([0 250]);
            xlim([ct(1) ct(end)]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)')
        end
    end
    sgtitle(sprintf('Direction: %d  Contrast: %0.2G', disp_direction, 1-disp_contrast));
end

%% within type Comparison
Fz = 100;
disp_direction = 180;
disp_contrast = 0;
disp_bar_witdth = [50, 100, 200, 400, 800];  % [50, 100, 200, 400, 800]
disp_speeds = [500, 1000, 2000, 4000, 8000]; % [500, 1000, 2000, 4000, 8000]
max_t = 459;
Disp_Type = 'ON';

cell_type_id = strcmpi(Disp_Type, 'ON');

ct = (0:max_t-1)/Fz;
Colors = parula(4);
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), location);
for i = 1:length(disp_contrast)
    Trace = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
    figure;
    for j = 1:length(disp_bar_witdth)
        for q = 1:length(disp_speeds)
            subplot(length(disp_bar_witdth), length(disp_speeds), (j-1)*length(disp_speeds)+q); hold on

            for k = 1:num_set
                dir_id = find(Data{k}.dim1_moving_direction == disp_direction);
                ctr_id = find(Data{k}.dim2_contrast == disp_contrast);
                bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j));
                sp_id = find(Data{k}.dim4_speeds == disp_speeds(q));
                csig = squeeze(mean(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
                Trace(j, q, k, :) = csig(1:max_t);
                if strcmpi(cell_type{k}, Disp_Type)
                    switch lower(location{k})
                        case 'temporal'
                            plot(ct, csig(1:max_t), 'Color', [180 0 180]/255);
                        case 'nasal'
                            plot(ct, csig(1:max_t), 'Color', [120 0 120]/255);
                    end
                end
            end
            gids_1 = cell_type_numeric==cell_type_id & location_type_numeric == 1;
            gids_2 = cell_type_numeric==cell_type_id & location_type_numeric == 0;
            plot(ct, squeeze(mean(Trace(j, q, gids_1, :), 3)), 'Color', [180 0 180]/255, 'LineWidth', 2)
            plot(ct, squeeze(mean(Trace(j, q, gids_2, :), 3)), 'Color', [120 0 120]/255, 'LineWidth', 2)
            ylim([0 250]);
            xlim([ct(1) ct(end)]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)')
            title(sprintf('%s group 1: n = %d, group 2: n = %d', Disp_Type, sum(gids_1), sum(gids_2)))
        end
    end
    sgtitle(sprintf('%s Direction: %d  Contrast: %0.2G', Disp_Type, disp_direction, 1-disp_contrast));
end

%% within type-location Comparison (4 figures: ON-Temporal, ON-Nasal, OFF-Temporal, OFF-Nasal)
save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';
Fz = 100;
disp_direction = 180;
disp_contrast = 0;
disp_bar_witdth = [100];  % [50, 100, 200, 400, 800]
disp_speeds = [500, 8000]; % [500, 1000, 2000, 4000, 8000]
max_t = 459;

ct = (0:max_t-1)/Fz;
group_colors = cat(3, ...
    [180, 0, 180; 120, 0, 120], ... % ON: Temporal, Nasal
    [0, 180, 0; 0, 120, 0]);        % OFF: Temporal, Nasal
group_colors = group_colors / 255;

cell_type_labels = {'ON', 'OFF'};
location_labels = {'Temporal', 'Nasal'};

cell_type_numeric = cellfun(@(x) strcmpi(x, 'ON'), cell_type);
location_type_numeric = cellfun(@(x) strcmpi(x, 'Temporal'), location);

for typeIdx = 1:2
    for locIdx = 1:2
        % Define current group
        Disp_Type = cell_type_labels{typeIdx};
        Disp_Location = location_labels{locIdx};
        % Logical indices for this group
        gids = (cell_type_numeric == (typeIdx==1)) & (location_type_numeric == (locIdx==1));
        
        if sum(gids)==0
            continue; % Skip if no cells in this group
        end

        color_this_group = squeeze(group_colors(locIdx, :, typeIdx));

        for i = 1:length(disp_contrast)
            Trace = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
            figure('Name', sprintf('%s-%s', Disp_Type, Disp_Location)); % Separate figure for each group
            
            for j = 1:length(disp_bar_witdth)
                for q = 1:length(disp_speeds)
                    subplot(length(disp_bar_witdth), length(disp_speeds), (j-1)*length(disp_speeds)+q); hold on
                    gids_idx = find(gids);
                    for idx = 1:length(gids_idx)
                        k = gids_idx(idx);
                        dir_id = find(Data{k}.dim1_moving_direction == disp_direction);
                        ctr_id = find(Data{k}.dim2_contrast == disp_contrast);
                        bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j));
                        sp_id = find(Data{k}.dim4_speeds == disp_speeds(q));
                        csig = squeeze(mean(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
                        Trace(j, q, k, :) = csig(1:max_t);
                        plot(ct, csig(1:max_t), 'Color', 0.5*ones(1, 3));
                    end
                    % Plot group mean
                    plot(ct, squeeze(mean(Trace(j, q, gids, :), 3)), 'Color', color_this_group, 'LineWidth', 2)
                    ylim([0 250]);
                    yticks(0:50:200);
                    yticklabels({'0', '', '100', '', '200'});
                    switch disp_speeds(q)
                        case 500
                            xlim([ct(1) 3.18]);
                            xticks(0:1:3);  
                            xticklabels({'0',  '1', '2',  '3'});
                        case 8000
                            xlim([ct(1) 1.13]);
                            xticks(0:0.4:1.2);  
                            xticklabels({'0',  '0.4', '0.8',  '1.2'});
                    end
                    
                    xlabel('Time (s)');
                    ylabel('Firing rate (spike/s)')
                    title(sprintf('%s-%s n = %d', Disp_Type, Disp_Location, sum(gids)))
                end
            end
            sgtitle(sprintf('%s-%s Direction: %d  Contrast: %0.2G', Disp_Type, Disp_Location, disp_direction, 1-disp_contrast));
            save_file_name = fullfile(save_folder, sprintf('MovingBar_%s-%s_%d_%0.2G', Disp_Type, Disp_Location, disp_direction, 1-disp_contrast));
            print(gcf, save_file_name, '-depsc', '-vector');
            print(gcf, save_file_name, '-dpng', '-r300');

        end
    end
end

function val = getfield_safe(s, fname)
    if isfield(s, fname)
        val = s.(fname);
    else
        val = nan;
    end
end