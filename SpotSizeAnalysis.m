% Analyze and plot size tuning (mean during stimulus ± SEM across cells)
% Requires: struct 'a' in workspace with fields:
%   - x1 (T x 1), y_labels (1 x S), and group arrays (T x (S * Ncells))
% Example groups: 'DN_ONSus_RF_UV', 'DN_ONSus_RF_GRN'
clc; close all;
save_fig_folder = './Figures/SpotSizeFit/';
if ~exist(save_fig_folder, 'dir')
    mkdir(save_fig_folder);
end
% ---------------- config ----------------
% Choose fitting method: 'fit_DoG_half_CDF' (default) or 'fitspotsizedog'
fitMethod = 'fitspotsizedog';  % Switch between 'fit_DoG_half_CDF' or 'fitspotsizedog'
% Choose model dimension for DoG half-CDF fitting: '1d' (default) or '2d'.
% Note: '2d' here means 2D line-accumulation (along a diameter), not disk area.
modelDim = '2d';  % set to '2d' to use 2D line-accumulation Gaussian half-CDF
if ~exist('a','var')
    error('Struct ''a'' not found in workspace.');
end
groups = {'AcuteZoneDT_ONSus_RF_GRN','DN_ONSus_RF_GRN'};  % edit/extend as needed
stim_idx = 101:300;  % during-stimulus period
sizes = a.y_labels(:)';                  % diameter labels
nSizes = numel(sizes);
T = numel(a.x1);

% Ensure sizes are ascending (plot from small to large)
[sorted_sizes, sort_idx] = sort(sizes, 'ascend');

% ---------------- compute ----------------
results = struct();
figure('Color','w'); hold on;
colororder = lines(numel(groups));

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

    % Mean across time during stimulus for each column (cell x size)
    stim_mean_by_col = mean(M(stim_idx, :), 1, 'omitnan');  % 1 x (nSizes*nCells)

    % Reshape into [nSizes x nCells]; columns order is per cell then sizes
    S = reshape(stim_mean_by_col, [nSizes, nCells]);

    % Reorder sizes to be ascending (and reorder S accordingly)
    S = S(sort_idx, :);

    % Mean and SEM across cells for each size
    mu = mean(S, 2, 'omitnan');                   % nSizes x 1
    nEff = sum(~isnan(S), 2);                     % per-size effective N
    sem = std(S, 0, 2, 'omitnan') ./ max(sqrt(nEff), 1);  % nSizes x 1

    % Store results
    results.(gname).sizes = sorted_sizes(:);
    results.(gname).mu = mu(:);
    results.(gname).sem = sem(:);
    results.(gname).perCell = S;  % size x cell matrix of per-cell means

    % Plot mean ± SEM
    h = errorbar(sorted_sizes, mu, sem, 'o-', ...
        'LineWidth', 1.5, 'Color', colororder(g,:), 'MarkerFaceColor', colororder(g,:));
    h.CapSize = 6;
    xlim([0, max(sorted_sizes)*1.05]);

    % ---------- Fit DoG model to mean curve ----------
    if strcmpi(fitMethod, 'fit_DoG_half_CDF')
        % Use original fit_DoG_half_CDF function
        [s_fit, y_fit, p_struct, fit_stats] = fit_DoG_half_CDF(sorted_sizes(:), mu(:), modelDim);
        results.(gname).fit.s_fit = s_fit;
        results.(gname).fit.y_fit = y_fit;
        results.(gname).fit.params = p_struct;
        results.(gname).fit.stats = fit_stats;
        
        % Display fit parameters for mean
        fprintf('\n%s DoG-CDF Fit Results (Mean, %s):\n', gname, p_struct.modelDim);
        fprintf('  Shared mean: mu=%.2f\n', p_struct.mu);
        fprintf('  Center sigma=%.2f, Surround sigma=%.2f\n', p_struct.sigma_c, p_struct.sigma_s);
        fprintf('  Model: gain=%.2f, w_s=%.3f, global_bias=%.2f\n', p_struct.gain, p_struct.w_s, p_struct.global_bias);
        fprintf('  Fit quality: loss=%.4f, R²=%.3f\n', fit_stats.fval, fit_stats.r_squared);
        
    elseif strcmpi(fitMethod, 'fitspotsizedog')
        % Use fitspotsizedog function
        fit_result = fitspotsizedog(sorted_sizes(:), mu(:), 'UseRobust', true);
        
        % Create dense curve for plotting
        s_fit = linspace(min(sorted_sizes), max(sorted_sizes), 100)';
        y_fit = fit_result.fun(s_fit);
        
        results.(gname).fit.s_fit = s_fit;
        results.(gname).fit.y_fit = y_fit;
        results.(gname).fit.params = fit_result.params;
        results.(gname).fit.stats = fit_result.gof;
        results.(gname).fit.center_weight_int = fit_result.center_weight_int;
        results.(gname).fit.surround_weight_int = fit_result.surround_weight_int;
        results.(gname).fit.weight_ratio = fit_result.weight_ratio;
        
        % Display fit parameters for mean
        fprintf('\n%s fitspotsizedog Fit Results (Mean):\n', gname);
        fprintf('  Center: k_c=%.2f, sigma_c=%.2f\n', fit_result.params(1), fit_result.params(2));
        fprintf('  Surround: k_s=%.2f, sigma_s=%.2f\n', fit_result.params(3), fit_result.params(4));
        fprintf('  Baseline: b=%.2f\n', fit_result.params(5));
        fprintf('  Weight ratio (surr/cent): %.3f\n', fit_result.weight_ratio);
        fprintf('  Fit quality: RMSE=%.4f, R²=%.3f\n', fit_result.gof.RMSE, fit_result.gof.R2);
    else
        error('Invalid fitMethod: %s. Use ''fit_DoG_half_CDF'' or ''fitspotsizedog''', fitMethod);
    end

    % Overlay fitted curve
    plot(s_fit, y_fit, '-', 'LineWidth', 2.0, 'Color', colororder(g,:));

    % ---------- Fit DoG model to each individual cell ----------
    fprintf('\nFitting individual cells for group %s...\n', gname);
    
    % Initialize storage for individual cell fits
    results.(gname).individual_fits = struct();
    cell_params = struct();
    
    % Set parameter fields based on fitting method
    if strcmpi(fitMethod, 'fit_DoG_half_CDF')
        param_fields = {'mu', 'sigma_c', 'sigma_s', 'w_s', 'gain', 'global_bias'};
        min_points = 5;  % Need at least 5 points for 6 parameters
    elseif strcmpi(fitMethod, 'fitspotsizedog')
        param_fields = {'k_c', 'sigma_c', 'k_s', 'sigma_s', 'baseline'};
        min_points = 4;  % Need at least 4 points for 5 parameters
    end
    
    % Configuration for spline interpolation and jittering
    n_interp_points = 10;  % Number of interpolated points between min and max size
    n_samples = 1;        % Number of jittered samples to generate
    jitter_std = 0.0;     % Standard deviation of jitter (relative to response range)
    
    for pf = 1:length(param_fields)
        cell_params.(param_fields{pf}) = nan(nCells, 1);
    end
    cell_r_squared = nan(nCells, 1);
    cell_rmse = nan(nCells, 1);
    
    for c = 1:nCells
        cell_response = S(:, c);  % Response of cell c to different sizes
        
        % Check if cell has enough valid data points
        valid_pts = ~isnan(cell_response);
        if sum(valid_pts) < min_points
            fprintf('  Cell %d: Insufficient data points (%d), skipping\n', c, sum(valid_pts));
            continue;
        end
        
        % Extract valid data for spline interpolation
        valid_sizes = sorted_sizes(valid_pts);
        valid_responses = cell_response(valid_pts);
        
        % Create spline interpolation and generate multiple jittered samples
        try
            % Create spline interpolation
            if length(valid_sizes) < 3
                % Use linear interpolation if less than 3 points
                interp_sizes = linspace(min(valid_sizes), max(valid_sizes), n_interp_points)';
                interp_responses = interp1(valid_sizes, valid_responses, interp_sizes, 'linear', 'extrap');
            else
                % Use spline interpolation
                interp_sizes = linspace(min(valid_sizes), max(valid_sizes), n_interp_points)';
                interp_responses = interp1(valid_sizes, valid_responses, interp_sizes, 'spline', 'extrap');
            end
            
            % Calculate jitter magnitude based on response range
            response_range = max(valid_responses) - min(valid_responses);
            jitter_magnitude = response_range * jitter_std;
            
            % Generate multiple jittered datasets
            all_sizes = [];
            all_responses = [];
            
            for sample = 1:n_samples
                % Add jitter to interpolated responses
                jittered_responses = interp_responses + jitter_magnitude * randn(size(interp_responses));
                
                % Accumulate all sampled data
                all_sizes = [all_sizes; interp_sizes];
                all_responses = [all_responses; jittered_responses];
            end
            
            % Fit using all accumulated data points
            if strcmpi(fitMethod, 'fit_DoG_half_CDF')
                % Fit DoG-CDF to aggregated spline-interpolated and jittered data
                [s_fit_cell, y_fit_cell, p_struct_cell, fit_stats_cell] = ...
                    fit_DoG_half_CDF(all_sizes, all_responses, modelDim);
                
                % Store individual cell fit results
                results.(gname).individual_fits.(sprintf('cell_%d', c)).s_fit = s_fit_cell;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).y_fit = y_fit_cell;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).params = p_struct_cell;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).stats = fit_stats_cell;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).original_sizes = valid_sizes;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).original_responses = valid_responses;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).n_samples_used = n_samples;
                
                % Store parameters for summary statistics
                for pf = 1:length(param_fields)
                    cell_params.(param_fields{pf})(c) = p_struct_cell.(param_fields{pf});
                end
                cell_r_squared(c) = fit_stats_cell.r_squared;
                cell_rmse(c) = fit_stats_cell.rmse;
                
                % Create parameter text for plot
                param_text = sprintf(['Spline+Jitter Fit (%dx samples):\n' ...
                    'Shared μ=%.2f\n' ...
                    'Center σ=%.2f, Surround σ=%.2f\n' ...
                    'Weight: w_s=%.3f\n' ...
                    'Gain: %.2f, Bias: %.2f'], ...
                    n_samples, p_struct_cell.mu, p_struct_cell.sigma_c, ...
                    p_struct_cell.sigma_s, ...
                    p_struct_cell.w_s, p_struct_cell.gain, p_struct_cell.global_bias);
                fit_title = sprintf('%s - Cell %d (%s, R²=%.3f, Spline+Jitter)', strrep(gname, '_', '\_'), c, p_struct_cell.modelDim, fit_stats_cell.r_squared);
                
            elseif strcmpi(fitMethod, 'fitspotsizedog')
                % Fit fitspotsizedog to aggregated spline-interpolated and jittered data
                fit_result_cell = fitspotsizedog(all_sizes, all_responses);
                
                % Create dense curve for plotting
                s_fit_cell = linspace(min(valid_sizes), max(valid_sizes), 100)';
                y_fit_cell = fit_result_cell.fun(s_fit_cell);
                
                % Store individual cell fit results
                results.(gname).individual_fits.(sprintf('cell_%d', c)).s_fit = s_fit_cell;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).y_fit = y_fit_cell;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).params = fit_result_cell.params;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).stats = fit_result_cell.gof;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).center_weight_int = fit_result_cell.center_weight_int;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).surround_weight_int = fit_result_cell.surround_weight_int;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).weight_ratio = fit_result_cell.weight_ratio;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).original_sizes = valid_sizes;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).original_responses = valid_responses;
                results.(gname).individual_fits.(sprintf('cell_%d', c)).n_samples_used = n_samples;
                
                % Store parameters for summary statistics
                cell_params.(param_fields{1})(c) = fit_result_cell.params(1); % k_c
                cell_params.(param_fields{2})(c) = fit_result_cell.params(2); % sigma_c
                cell_params.(param_fields{3})(c) = fit_result_cell.params(3); % k_s
                cell_params.(param_fields{4})(c) = fit_result_cell.params(4); % sigma_s
                cell_params.(param_fields{5})(c) = fit_result_cell.params(5); % baseline
                cell_r_squared(c) = fit_result_cell.gof.R2;
                cell_rmse(c) = fit_result_cell.gof.RMSE;
                
                % Create parameter text for plot
                param_text = sprintf(['Spline+Jitter Fit (%dx samples):\n' ...
                    'Center: k_c=%.2f, σ_c=%.2f\n' ...
                    'Surround: k_s=%.2f, σ_s=%.2f\n' ...
                    'Baseline: %.2f\n' ...
                    'Weight ratio: %.3f'], ...
                    n_samples, fit_result_cell.params(1), fit_result_cell.params(2), ...
                    fit_result_cell.params(3), fit_result_cell.params(4), ...
                    fit_result_cell.params(5), fit_result_cell.weight_ratio);
                fit_title = sprintf('%s - Cell %d (fitspotsizedog, R²=%.3f, Spline+Jitter)', strrep(gname, '_', '\_'), c, fit_result_cell.gof.R2);
            end
            
            % Create and save individual cell plot
            fig_cell = figure('Color', 'w', 'Position', [100, 100, 800, 500]);
            subplot(1, 2, 1);
            hold on;
            
            % Plot original data points
            plot(valid_sizes, valid_responses, 'o', 'MarkerSize', 10, ...
                'Color', colororder(g,:), 'MarkerFaceColor', colororder(g,:), ...
                'LineWidth', 2, 'DisplayName', 'Original Data');
            % Plot fitted curve
            plot(s_fit_cell, y_fit_cell, '-', 'LineWidth', 3, 'Color', [0.8 0.2 0.2], ...
                'DisplayName', 'Fitted Curve');
            % Styling
            grid on;
            xlim([0, max(sorted_sizes)*1.05]);
            xlabel('Spot Diameter');
            ylabel('Mean Firing Rate (Hz)');
            title('Original Data & Fit');
            legend('Location', 'best');
            
            % Second subplot: show some of the interpolated/jittered data
            subplot(1, 2, 2);
            hold on;
            % Show a subset of jittered samples (first 3 to avoid clutter)
            for sample = 1:min(3, n_samples)
                start_idx = (sample-1)*n_interp_points + 1;
                end_idx = sample*n_interp_points;
                plot(all_sizes(start_idx:end_idx), all_responses(start_idx:end_idx), '.', ...
                    'Color', [0.7 0.7 0.7], 'MarkerSize', 4, 'HandleVisibility', 'off');
            end
            % Plot original data points on top
            plot(valid_sizes, valid_responses, 'o', 'MarkerSize', 10, ...
                'Color', colororder(g,:), 'MarkerFaceColor', colororder(g,:), ...
                'LineWidth', 2, 'DisplayName', 'Original Data');
            % Plot fitted curve
            plot(s_fit_cell, y_fit_cell, '-', 'LineWidth', 3, 'Color', [0.8 0.2 0.2], ...
                'DisplayName', 'Fitted Curve');
            grid on;
            xlim([0, max(sorted_sizes)*1.05]);
            xlabel('Spot Diameter');
            ylabel('Mean Firing Rate (Hz)');
            title('Spline+Jitter Samples');
            legend('Location', 'best');
            
            % Add overall title and parameter text
            sgtitle(fit_title);
            % Add parameter text to the figure
            annotation('textbox', [0.80, 0.20, 0.3, 0.15], 'String', param_text, ...
                'FontSize', 9, 'BackgroundColor', [1 1 1 0.9], 'EdgeColor', 'k', ...
                'FitBoxToText', 'on');
            
            % Save figure
            filename = sprintf('%s_Cell_%d_SpotSizeFit_SplineJitter.png', gname, c);
            saveas(fig_cell, fullfile(save_fig_folder, filename));
            close(fig_cell);
            
            if strcmpi(fitMethod, 'fit_DoG_half_CDF')
                fprintf('  Cell %d: R²=%.3f, RMSE=%.3f (using %d interpolated samples), saved as %s\n', ...
                    c, fit_stats_cell.r_squared, fit_stats_cell.rmse, n_samples*n_interp_points, filename);
            elseif strcmpi(fitMethod, 'fitspotsizedog')
                fprintf('  Cell %d: R²=%.3f, RMSE=%.3f (using %d interpolated samples), saved as %s\n', ...
                    c, fit_result_cell.gof.R2, fit_result_cell.gof.RMSE, n_samples*n_interp_points, filename);
            end
            
        catch ME
            fprintf('  Cell %d: Fit failed - %s\n', c, ME.message);
        end
    end
    
    % Store summary statistics of individual cell parameters
    results.(gname).cell_params_summary = struct();
    for pf = 1:length(param_fields)
        valid_params = ~isnan(cell_params.(param_fields{pf}));
        if any(valid_params)
            results.(gname).cell_params_summary.(param_fields{pf}).mean = ...
                mean(cell_params.(param_fields{pf})(valid_params));
            results.(gname).cell_params_summary.(param_fields{pf}).std = ...
                std(cell_params.(param_fields{pf})(valid_params));
            results.(gname).cell_params_summary.(param_fields{pf}).values = ...
                cell_params.(param_fields{pf});
        end
    end
    
    % Store fit quality summary
    valid_r2 = ~isnan(cell_r_squared);
    if any(valid_r2)
        results.(gname).fit_quality.r_squared_mean = mean(cell_r_squared(valid_r2));
        results.(gname).fit_quality.r_squared_std = std(cell_r_squared(valid_r2));
        results.(gname).fit_quality.rmse_mean = mean(cell_rmse(valid_r2));
        results.(gname).fit_quality.rmse_std = std(cell_rmse(valid_r2));
        results.(gname).fit_quality.n_successful_fits = sum(valid_r2);
    end
    
    fprintf('Completed individual cell fitting for %s: %d/%d cells successfully fit\n', ...
        gname, sum(valid_r2), nCells);
end

% ---------------- figure styling ----------------
grid on;
xlim([0, max(sorted_sizes)*1.05]);
xlabel('Diameter');
ylabel('Mean firing rate during stimulus (Hz)');
if strcmpi(fitMethod, 'fit_DoG_half_CDF')
    title(sprintf('Size tuning: mean ± SEM across cells with DoG-CDF fit (%s)', upper(modelDim)));
elseif strcmpi(fitMethod, 'fitspotsizedog')
    title('Size tuning: mean ± SEM across cells with fitspotsizedog fit');
end
legend(groups, 'Interpreter', 'none', 'Location', 'best');
hold off;

% Save the summary figure
if strcmpi(fitMethod, 'fit_DoG_half_CDF')
    suffix = upper(modelDim);
    summary_filename = sprintf('SpotSizeAnalysis_Summary_%s.png', suffix);
elseif strcmpi(fitMethod, 'fitspotsizedog')
    summary_filename = 'SpotSizeAnalysis_Summary_fitspotsizedog.png';
end
saveas(gcf, fullfile(save_fig_folder, summary_filename));
fprintf('\nSummary figure saved as: %s\n', summary_filename);

% Print summary statistics for all groups
fprintf('\n========== SUMMARY STATISTICS ==========\n');
for g = 1:numel(groups)
    gname = groups{g};
    if isfield(results, gname) && isfield(results.(gname), 'fit_quality')
        fprintf('\n%s:\n', gname);
        fprintf('  Successfully fit: %d cells\n', results.(gname).fit_quality.n_successful_fits);
        fprintf('  Mean R²: %.3f ± %.3f\n', ...
            results.(gname).fit_quality.r_squared_mean, results.(gname).fit_quality.r_squared_std);
        fprintf('  Mean RMSE: %.3f ± %.3f\n', ...
            results.(gname).fit_quality.rmse_mean, results.(gname).fit_quality.rmse_std);
        
        % Print parameter statistics
        if isfield(results.(gname), 'cell_params_summary')
            if strcmpi(fitMethod, 'fit_DoG_half_CDF')
                param_fields_print = {'mu', 'sigma_c', 'sigma_s', 'w_s', 'gain', 'global_bias'};
                param_names = {'Shared μ', 'Center σ', 'Surround σ', 'Surround weight', 'Gain', 'Bias'};
            elseif strcmpi(fitMethod, 'fitspotsizedog')
                param_fields_print = {'k_c', 'sigma_c', 'k_s', 'sigma_s', 'baseline'};
                param_names = {'Center gain', 'Center σ', 'Surround gain', 'Surround σ', 'Baseline'};
            end
            
            for pf = 1:length(param_fields_print)
                if isfield(results.(gname).cell_params_summary, param_fields_print{pf})
                    param_stats = results.(gname).cell_params_summary.(param_fields_print{pf});
                    fprintf('  %s: %.3f ± %.3f\n', param_names{pf}, param_stats.mean, param_stats.std);
                end
            end
        end
    end
end


%%
close all; clc;
% Step 1: Define the x-range
x = 0:0.1:5;

% Define mean and standard deviation for the normal distribution
mu = 0; % mean
sigma = 1; % standard deviation

% Step 2: Calculate the CDF values
y = normcdf(x, mu, sigma);
y2 = normcdf(x, mu, sigma*1.5);
figure; hold on
% Step 3: Plot the CDF
plot(x, y);
plot(x, y-0.5*y2);
title('Theoretical Normal CDF');
xlabel('x');
ylabel('F(x)');
ylim([0 1]);
grid on;