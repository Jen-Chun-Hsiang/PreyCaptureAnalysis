close all; clear; clc;
%% Load processed data
set_name = 'latest';  % before081425
is_show_fitted = 1;

%% ============ USER PARAMETERS ============
Disp_Type = 'ON';  % 'ON' or 'OFF' - cell type to display in "within type Comparison" 5x5 grid figure
is_show_rf_edge_timeline = true;  % Add dashed vertical lines for bar-edge entry/exit at RF center
rf_center_diameter_um = 300;
moving_bar_startposition_um = 500; % IN.startposition in StageVSS script (um)

% d101924 unlikely to be OFF temporal, responses is very quick to sustained parts

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
                    'b101924', 'c101924',            'e101924', 'b103124',   'e103124', 'a110424',...
                    'c110424',                       'f110424', 'g110424',   'a110924', 'b110924', 'c110924',...
                    'a111224'}; % removed d101924 'OFF' 'Nasal'
        cell_type = {'OFF',      'OFF',    'OFF',      'ON',       'OFF',     'ON',      'OFF',...
                    'OFF',      'OFF',    'ON',       'ON',       'ON',      'ON',      'ON',       'ON',...
                    'ON',       'OFF',                'OFF',      'ON',      'OFF',     'ON',...
                    'ON',                             'ON',       'OFF',     'ON',      'OFF',      'OFF',...
                    'ON'};
        location =  {'Temporal', 'Temporal','Nasal',   'Nasal',    'Nasal',   'Nasal',   'Nasal',...
                    'Temporal', 'Temporal','Temporal','Temporal', 'Nasal',   'Nasal',   'Nasal',    'Nasal',...
                    'Nasal',    'Nasal',             'Nasal',    'Temporal','Temporal','Temporal', ...
                    'Temporal',                      'Temporal',  'Temporal','Temporal','Temporal', 'Temporal',...
                    'Temporal'};
end

save_fig_folder = './Figures/MovingBarSummary/';
if ~exist(save_fig_folder, 'dir')
    mkdir(save_fig_folder);
end
clear Data FitD
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
         'BaselineCorr', 'LNK_params_s', 'LNK_params_w', 'LN_params_s', 'LNK_params_d', 'sim_nl_s', 'trail', 'trial_table');
        Cdat(i, :) = loaded_data.PredictionResults;
        Cbas(i, :) = loaded_data.BaselineCorr;
        Csw(i, :) = [loaded_data.LNK_params_s.w_xs loaded_data.LNK_params_w.w_xs,...
                     loaded_data.LN_params_s.gamma loaded_data.LNK_params_d.w_xs];
        FitD{i}.sim_nl_s = loaded_data.sim_nl_s;
        FitD{i}.trail = loaded_data.trail;
        FitD{i}.trial_table = loaded_data.trial_table;

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
%% check FitD.trail and FitD.trial_table consistency
check_ids = randperm(num_set, 3)
% assert(numel(FitD{check_ids(1)}.trail) == numel(FitD{check_ids(2)}.trail));
assert(sum(FitD{check_ids(1)}.trial_table - FitD{check_ids(3)}.trial_table, 'all') == 0);
% assert(sum(FitD{check_ids(2)}.trial - FitD{check_ids(3)}.trial, 'all') == 0);


%%
keyboard;


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

    % ---- Nonparametric paired comparison: Repeat reliability vs LN+surround ----
    % Each data point is a cell within each group.
    % Uses Wilcoxon signed-rank test (signrank) on paired samples.
    np_results = struct();
    for g = 1:4
        group_mask = group_indices{g};
        base = Cbas(group_mask);
        lns = Cdat(group_mask, 9); % 9th column plotted as LN_s in this script

        valid = ~isnan(base) & ~isnan(lns);
        base = base(valid);
        lns = lns(valid);

        res = struct();
        res.group = groups{g};
        res.n = numel(base);
        res.median_repeat_reliability = median(base, 'omitnan');
        res.median_LN_surround = median(lns, 'omitnan');
        res.median_diff_LNsurround_minus_repeat = median(lns - base, 'omitnan');

        if res.n >= 3
            try
                [p, h, stats] = signrank(lns, base); % two-sided by default
                res.test = 'signrank (paired, two-sided)';
                res.p = p;
                res.h = h;
                res.signedrank = getfield_safe(stats, 'signedrank');
                res.zval = getfield_safe(stats, 'zval');
            catch ME
                res.test = 'signrank failed';
                res.p = nan;
                res.h = nan;
                res.error = ME.message;
            end
        else
            res.test = 'signrank skipped (n<3)';
            res.p = nan;
            res.h = nan;
        end

        np_results.(matlab.lang.makeValidName(groups{g})) = res;
    end

    fprintf('\n=== NONPARAMETRIC PAIRED COMPARISON: Repeat reliability vs LN+surround ===\n');
    for g = 1:4
        key = matlab.lang.makeValidName(groups{g});
        res = np_results.(key);
        fprintf('%s: n=%d, median(RR)=%.3f, median(LN_s)=%.3f, median(diff)=%.3f, p=%g\n', ...
            res.group, res.n, res.median_repeat_reliability, res.median_LN_surround, res.median_diff_LNsurround_minus_repeat, res.p);
    end

    % Save stats next to other summary outputs
    try
        save(fullfile(save_fig_folder, 'NonparametricComparison_RR_vs_LNs.mat'), 'np_results');
    catch
        % If save_fig_folder is not writable for some reason, fall back to CWD
        save('NonparametricComparison_RR_vs_LNs.mat', 'np_results');
    end
    
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
        plot(x_vals, sorted_cdat_lnk, 'r-', 'LineWidth', 2, 'DisplayName', 'LN_s');
        
        xlabel('Cell (sorted by repeat reliability)');
        ylabel('Correlation coefficient');
        title([groups{g} ' (n=' num2str(sum(group_mask)) ')']);
        xticks(1:2:length(sorted_cbas));
        xticklabels(arrayfun(@num2str, 1:2:length(sorted_cbas), 'UniformOutput', false));
        yticks(0:0.25:1);
        yticklabels({'0', '', '0.5', '', '1'});
        ylim([0 1]);
        xlim([1 length(sorted_cbas)]);
        xticks(1:2:length(sorted_cbas));
        xticklabels(arrayfun(@num2str, 1:2:length(sorted_cbas), 'UniformOutput', false));
        legend('Location', 'best');
        grid off;
    end
    sgtitle('Model Performance: LN vs LNK vs Repeat Reliability');

    saveas(gcf, fullfile(save_fig_folder, 'Model-Performance-LN-LNK.png'));
    print(gcf, fullfile(save_fig_folder, 'Model-Performance-LN-LNK'), '-depsc', '-vector');

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
disp_direction = 0;
disp_contrast = 0;
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
    sgtitle(sprintf('ON vs OFF Comparison - Grid: Bar Width (rows) × Speed (cols)\nDirection: %d  Contrast: %0.2G', disp_direction, 1-disp_contrast));
end

%% within type Comparison
Fz = 100;
disp_direction = 0;
disp_contrast = 0;
disp_bar_witdth = [50, 100, 200, 400, 800];  % [50, 100, 200, 400, 800]
disp_speeds = [500, 1000, 2000, 4000, 8000]; % [500, 1000, 2000, 4000, 8000]
max_t = 459;

cell_type_id = strcmpi(Disp_Type, 'ON');

ct = (0:max_t-1)/Fz;
Colors = parula(4);
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), location);

if strcmpi(Disp_Type, 'ON')
    temporal_mean_color = [180 0 180]/255;
    nasal_mean_color = [120 0 120]/255;
else
    temporal_mean_color = [0 180 0]/255;
    nasal_mean_color = [0 120 0]/255;
end
lighten_factor = 0.5;
temporal_ind_color = temporal_mean_color + lighten_factor*(1 - temporal_mean_color);
nasal_ind_color = nasal_mean_color + lighten_factor*(1 - nasal_mean_color);

rf_radius_um = rf_center_diameter_um/2;
save_folder_fig2 = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';
if ~exist(save_folder_fig2, 'dir')
    mkdir(save_folder_fig2);
end

for i = 1:length(disp_contrast)
    Trace = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
    figure('Color', 'w', 'Position', [100, 80, 1350, 950]);
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
                            plot(ct, csig(1:max_t), 'Color', temporal_ind_color, 'LineWidth', 0.8);
                        case 'nasal'
                            plot(ct, csig(1:max_t), 'Color', nasal_ind_color, 'LineWidth', 0.8);
                    end
                end
            end
            gids_1 = cell_type_numeric==cell_type_id & location_type_numeric == 1;
            gids_2 = cell_type_numeric==cell_type_id & location_type_numeric == 0;
            plot(ct, squeeze(mean(Trace(j, q, gids_1, :), 3)), 'Color', temporal_mean_color, 'LineWidth', 2.5)
            plot(ct, squeeze(mean(Trace(j, q, gids_2, :), 3)), 'Color', nasal_mean_color, 'LineWidth', 2.5)

            if is_show_rf_edge_timeline
                bar_center_start_um = moving_bar_startposition_um + 0.5*disp_bar_witdth(j);
                edge_offset_um = rf_radius_um + 0.5*disp_bar_witdth(j);
                t_edge_enter = (bar_center_start_um - edge_offset_um) / disp_speeds(q);
                t_edge_exit = (bar_center_start_um + edge_offset_um) / disp_speeds(q);
                xline(t_edge_enter, '--', 'Color', 0.2*[1 1 1], 'LineWidth', 1.1);
                xline(t_edge_exit, '--', 'Color', 0.2*[1 1 1], 'LineWidth', 1.1);
            end

            ylim([0 250]);
            xlim([0 4]);
            xticks([0 2 4]);
            xticklabels({'0', '2', '4'});
            yticks([0 100 200]);
            set(gca, 'Box', 'off', 'TickDir', 'out', 'LineWidth', 1, 'FontName', 'Arial', 'FontSize', 9, 'Layer', 'top');
            if j < length(disp_bar_witdth)
                xticklabels({});
            else
                xlabel('Time (s)', 'FontWeight', 'bold');
            end
            if q > 1
                yticklabels({});
            else
                ylabel('Firing rate (spike/s)', 'FontWeight', 'bold');
            end
            if j == 1
                title(sprintf('Speed %d', disp_speeds(q)), 'FontWeight', 'bold');
            end
            if q == 1
                text(-0.55, 125, sprintf('W=%d', disp_bar_witdth(j)), ...
                    'Rotation', 90, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'FontName', 'Arial', 'FontSize', 9, 'FontWeight', 'bold');
            end
            if j == 1 && q == 1 && is_show_rf_edge_timeline
                text(0.06, 232, sprintf('--: edge @ RF %d\x03BCm, IN.startposition=%d\x03BCm', rf_center_diameter_um, moving_bar_startposition_um), ...
                    'FontName', 'Arial', 'FontSize', 8, 'Color', [0.2 0.2 0.2]);
            end
        end
    end
    sgtitle(sprintf('%s Cells: Temporal (light) vs Nasal (dark) | Direction %d | Contrast %0.2G\nRows: bar width (\x03BCm), Columns: speed | Dashed lines: bar-edge entry/exit at RF center', Disp_Type, disp_direction, 1-disp_contrast), 'FontWeight', 'bold');
    save_file_name = fullfile(save_folder_fig2, sprintf('MovingBar_WithinType_5x5_%s_Dir%d_Contrast%0.2G', Disp_Type, disp_direction, 1-disp_contrast));
    print(gcf, save_file_name, '-depsc', '-painters'); % EPS format
    print(gcf, save_file_name, '-dpng', '-r300'); % PNG, 600 dpi
end

%% within type-location Comparison (4 figures: ON-Temporal, ON-Nasal, OFF-Temporal, OFF-Nasal)
save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';
Fz = 100;
disp_direction = 0;
disp_contrast = 0;
disp_bar_witdth = [800];  % [50, 100, 200, 400, 800]
disp_speeds = [1000 8000]; % [500, 1000, 2000, 4000, 8000]
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
            Trace_s = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
            figure('Name', sprintf('%s-%s', Disp_Type, Disp_Location)); % Separate figure for each group
            
            for j = 1:length(disp_bar_witdth)
                for q = 1:length(disp_speeds)
                    subplot(length(disp_speeds), length(disp_bar_witdth), (q-1)*length(disp_bar_witdth)+j); hold on
                    gids_idx = find(gids);
                    for idx = 1:length(gids_idx)
                        k = gids_idx(idx);
                        dir_id = find(Data{k}.dim1_moving_direction == disp_direction);
                        ctr_id = find(Data{k}.dim2_contrast == disp_contrast);
                        bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j));
                        sp_id = find(Data{k}.dim4_speeds == disp_speeds(q));
                        csig = squeeze(mean(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
                        csim_tab  = getfield_safe(FitD{k}, 'trial_table');
                        if ~isnan(csim_tab)
                            sim =  FitD{k}.sim_nl_s;
                            ctrail = FitD{k}.trail;
                            tids = find(csim_tab(:,1)==ctr_id & csim_tab(:,2)==bw_id & csim_tab(:,3)==sp_id);
                            clear csim_c csim
                            csim_c = cell(1, length(tids)); 
                            for  tt = 1:length(tids)
                                tid = tids(tt);
                                csim_c{tt} = sim(ctrail == tid);
                            end
                            temp_Ls = max(cellfun(@(x) length(x), csim_c));
                            csim = nan(length(tids), temp_Ls);
                            for tt = 1:length(tids)
                                csim(tt, 1:length(csim_c{tt})) = csim_c{tt};
                            end
                            csim = mean(csim, 1);
                            if length(csim) < max_t
                                csim = [csim, nan(1, max_t-length(csim))];
                            end
                            Trace_s(j, q, k, :) = csim(1:max_t);
                            
                        end

                        
                        Trace(j, q, k, :) = csig(1:max_t);
                        % plot(ct, csig(1:max_t), 'Color', 0.5*ones(1, 3));
                    end
                    % Plot group mean
                    mean_trace = squeeze(mean(Trace(j, q, gids, :), 3));
                    sem_trace = squeeze(std(Trace(j, q, gids, :), 0, 3)) / sqrt(sum(gids));
                    mean_trace_s = squeeze(mean(Trace_s(j, q, gids, :), 3));
                    sem_trace_s = squeeze(std(Trace_s(j, q, gids, :), 0, 3)) / sqrt(sum(gids));
                    % Plot shaded error bars for mean ± SEM
                    shadePlot(ct, mean_trace, sem_trace,  color_this_group);
                    shadePlot(ct, mean_trace_s, sem_trace_s, 'k');
                    ylim([0 250]);
                    yticks(0:50:200);
                    yticklabels({'0', '', '100', '', '200'});
                    switch disp_speeds(q)
                        case 500
                            xlim([ct(1) 3.18]);
                            xticks(0:1:3);  
                            xticklabels({'0',  '1', '2',  '3'});
                        case 1000
                            xlim([ct(1) 2.77]);
                            xticks(0:1:2);  
                            xticklabels({'0',  '1', '2'});
                        case 8000
                            xlim([ct(1) 1.2]);
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

%% Quantify peak responses across speed and bar width (group heatmaps)
Fz = 100;
disp_direction = 0;
disp_contrast = 0;
disp_bar_witdth = [50, 100, 200, 400, 800];
disp_speeds = [500, 1000, 2000, 4000, 8000];
max_t = 459;
is_baseline_substraction = false;
baseline_percentile = 0.05; % Use [0,1], e.g. 0.05=5th percentile, 0.5=median
save_folder_heatmap = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';
if ~exist(save_folder_heatmap, 'dir')
    mkdir(save_folder_heatmap);
end

if baseline_percentile < 0 || baseline_percentile > 1
    error('baseline_percentile must be between 0 and 1. Example: 0.05 for 5%%, 0.5 for median.');
end

if is_baseline_substraction
    metric_label = 'BaselineSubtractedPeak';
    metric_title = sprintf('peak - %.1f%% quantile', baseline_percentile * 100);
else
    metric_label = 'RawPeak';
    metric_title = 'raw peak';
end

cell_type_numeric = cellfun(@(x) strcmpi(x, 'ON'), cell_type);
location_type_numeric = cellfun(@(x) strcmpi(x, 'Temporal'), location);

group_names = {'ON-Temporal', 'ON-Nasal', 'OFF-Temporal', 'OFF-Nasal'};
group_masks = {
    cell_type_numeric == 1 & location_type_numeric == 1; ...
    cell_type_numeric == 1 & location_type_numeric == 0; ...
    cell_type_numeric == 0 & location_type_numeric == 1; ...
    cell_type_numeric == 0 & location_type_numeric == 0  ...
};

n_width = numel(disp_bar_witdth);
n_speed = numel(disp_speeds);

peak_by_cell = nan(n_width, n_speed, num_set);

for k = 1:num_set
    dir_id = find(Data{k}.dim1_moving_direction == disp_direction, 1);
    ctr_id = find(Data{k}.dim2_contrast == disp_contrast, 1);
    if isempty(dir_id) || isempty(ctr_id)
        continue;
    end

    for j = 1:n_width
        bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j), 1);
        if isempty(bw_id)
            continue;
        end

        for q = 1:n_speed
            sp_id = find(Data{k}.dim4_speeds == disp_speeds(q), 1);
            if isempty(sp_id)
                continue;
            end

            trial_data = squeeze(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :));
            peak_by_cell(j, q, k) = compute_recording_peak_metric(trial_data, max_t, is_baseline_substraction, baseline_percentile);
        end
    end
end

% Temporal vs Nasal statistics at each width-speed cell (within ON and within OFF)
% Type index: 1=ON, 2=OFF
p_raw_by_type = nan(n_width, n_speed, 2);
p_fdr_by_type = nan(n_width, n_speed, 2);
sig_markers_by_type = cell(n_width, n_speed, 2);
n_temporal_by_type = zeros(n_width, n_speed, 2);
n_nasal_by_type = zeros(n_width, n_speed, 2);

for type_idx = 1:2
    if type_idx == 1
        temporal_mask = cell_type_numeric == 1 & location_type_numeric == 1;
        nasal_mask = cell_type_numeric == 1 & location_type_numeric == 0;
    else
        temporal_mask = cell_type_numeric == 0 & location_type_numeric == 1;
        nasal_mask = cell_type_numeric == 0 & location_type_numeric == 0;
    end

    p_vec = nan(n_width*n_speed, 1);
    vec_idx = 0;
    for j = 1:n_width
        for q = 1:n_speed
            vec_idx = vec_idx + 1;
            vals_temporal = squeeze(peak_by_cell(j, q, temporal_mask));
            vals_nasal = squeeze(peak_by_cell(j, q, nasal_mask));
            vals_temporal = vals_temporal(~isnan(vals_temporal));
            vals_nasal = vals_nasal(~isnan(vals_nasal));
            n_temporal_by_type(j, q, type_idx) = numel(vals_temporal);
            n_nasal_by_type(j, q, type_idx) = numel(vals_nasal);

            if numel(vals_temporal) >= 2 && numel(vals_nasal) >= 2
                p_raw_by_type(j, q, type_idx) = ranksum(vals_temporal, vals_nasal);
                p_vec(vec_idx) = p_raw_by_type(j, q, type_idx);
            end
        end
    end

    p_fdr_vec = fdr_bh_builtin(p_vec);
    vec_idx = 0;
    for j = 1:n_width
        for q = 1:n_speed
            vec_idx = vec_idx + 1;
            p_fdr_by_type(j, q, type_idx) = p_fdr_vec(vec_idx);
            sig_markers_by_type{j, q, type_idx} = p_to_marker(p_fdr_vec(vec_idx));
        end
    end
end

% Console report for temporal vs nasal stats
fprintf('\n=== Temporal vs Nasal stats on heatmap cells (%s) ===\n', metric_label);
for type_idx = 1:2
    if type_idx == 1
        type_label = 'ON';
    else
        type_label = 'OFF';
    end
    fprintf('\n[%s] width_um  speed_um_s  n_temporal  n_nasal    p_raw      p_fdr   sig\n', type_label);
    for j = 1:n_width
        for q = 1:n_speed
            p_raw_val = p_raw_by_type(j, q, type_idx);
            p_fdr_val = p_fdr_by_type(j, q, type_idx);
            sig_mark = sig_markers_by_type{j, q, type_idx};
            if isempty(sig_mark)
                sig_mark = '';
            end
            fprintf('%8d %11d %11d %8d %10.3g %10.3g   %s\n', ...
                disp_bar_witdth(j), disp_speeds(q), ...
                n_temporal_by_type(j, q, type_idx), n_nasal_by_type(j, q, type_idx), ...
                p_raw_val, p_fdr_val, sig_mark);
        end
    end
end

peak_heatmap_group = nan(n_width, n_speed, numel(group_names));
peak_count_group = zeros(n_width, n_speed, numel(group_names));

for g = 1:numel(group_names)
    gids = group_masks{g};
    if ~any(gids)
        continue;
    end

    for j = 1:n_width
        for q = 1:n_speed
            vals = squeeze(peak_by_cell(j, q, gids));
            peak_heatmap_group(j, q, g) = mean(vals, 'omitnan');
            peak_count_group(j, q, g) = sum(~isnan(vals));
        end
    end

    fig = figure('Name', sprintf('Peak heatmap %s (%s)', group_names{g}, metric_label));
    imagesc(peak_heatmap_group(:, :, g));
    set(gca, 'YDir', 'normal');
    colormap(parula);
    cb = colorbar;
    caxis([50 210]);
    cb.Ticks = [50, 100, 150, 200];
    cb.TickLabels = {'50', '100', '150', '200'};

    xticks(1:n_speed);
    xticklabels(arrayfun(@num2str, disp_speeds, 'UniformOutput', false));
    yticks(1:n_width);
    yticklabels(arrayfun(@num2str, disp_bar_witdth, 'UniformOutput', false));

    xlabel('Speed');
    ylabel('Bar width');
    title(sprintf('%s: mean %s (n=%d)', group_names{g}, metric_title, sum(gids)));

    % Overlay significance markers (Temporal vs Nasal, FDR-corrected)
    if g <= 2
        type_idx = 1; % ON
    else
        type_idx = 2; % OFF
    end
    for j = 1:n_width
        for q = 1:n_speed
            mark = sig_markers_by_type{j, q, type_idx};
            if ~isempty(mark)
                text(q, j, mark, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
                    'Color', 'k', 'FontWeight', 'bold', 'FontSize', 11);
            end
        end
    end

    save_name = sprintf('MovingBar_PeakHeatmap_%s_%s', metric_label, strrep(group_names{g}, '-', '_'));
    save_file_path = fullfile(save_folder_heatmap, save_name);
    print(gcf, save_file_path, '-depsc', '-painters'); % EPS format
    print(gcf, save_file_path, '-dpng', '-r300'); % PNG, 300 dpi
end

save(fullfile(save_fig_folder, sprintf('MovingBar_PeakHeatmap_Data_%s.mat', metric_label)), ...
    'peak_by_cell', 'peak_heatmap_group', 'peak_count_group', ...
    'p_raw_by_type', 'p_fdr_by_type', 'sig_markers_by_type', ...
    'n_temporal_by_type', 'n_nasal_by_type', ...
    'disp_bar_witdth', 'disp_speeds', 'group_names', 'group_masks', ...
    'disp_direction', 'disp_contrast', 'max_t', 'Fz', ...
    'is_baseline_substraction', 'baseline_percentile', 'metric_label', 'metric_title');

function recording_metric = compute_recording_peak_metric(trial_data, max_t, is_baseline_substraction, baseline_percentile)
    if isempty(trial_data)
        recording_metric = nan;
        return;
    end

    if isvector(trial_data)
        traces = trial_data(:)';
    else
        if size(trial_data, 1) <= size(trial_data, 2)
            traces = trial_data;
        else
            traces = trial_data';
        end
    end

    trace_end = min(max_t, size(traces, 2));
    if trace_end < 1
        recording_metric = nan;
        return;
    end

    traces = traces(:, 1:trace_end);
    trial_metrics = nan(size(traces, 1), 1);

    for t = 1:size(traces, 1)
        tr = traces(t, :);
        tr = tr(~isnan(tr));
        if isempty(tr)
            continue;
        end

        peak_val = max(tr);
        if is_baseline_substraction
            baseline_val = prctile(tr, baseline_percentile * 100);
            trial_metrics(t) = peak_val - baseline_val;
        else
            trial_metrics(t) = peak_val;
        end
    end

    recording_metric = mean(trial_metrics, 'omitnan');
end

function p_fdr = fdr_bh_builtin(p_raw)
    p_fdr = nan(size(p_raw));
    valid = ~isnan(p_raw);
    if ~any(valid)
        return;
    end

    if exist('mafdr', 'file') == 2
        p_fdr(valid) = mafdr(p_raw(valid), 'BHFDR', true);
    else
        warning('mafdr not found. Falling back to internal BH implementation.');
        p_fdr = fdr_bh_fallback(p_raw);
    end
end

function p_fdr = fdr_bh_fallback(p_raw)
    p_fdr = nan(size(p_raw));
    valid = ~isnan(p_raw);
    if ~any(valid)
        return;
    end

    pv = p_raw(valid);
    m = numel(pv);
    [p_sorted, order] = sort(pv(:));
    q_sorted = p_sorted .* (m ./ (1:m)');
    q_sorted = min(1, q_sorted);

    for i = m-1:-1:1
        q_sorted(i) = min(q_sorted(i), q_sorted(i+1));
    end

    qv = nan(size(pv));
    qv(order) = q_sorted;
    p_fdr(valid) = qv;
end

function marker = p_to_marker(p_fdr)
    marker = '';
    if isnan(p_fdr)
        return;
    end

    if p_fdr < 0.001
        marker = '***';
    elseif p_fdr < 0.01
        marker = '**';
    elseif p_fdr < 0.05
        marker = '*';
    end
end

function val = getfield_safe(s, fname)
    if isfield(s, fname)
        val = s.(fname);
    else
        val = nan;
    end
end