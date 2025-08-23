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
if ~exist('a','var')
    error('Struct ''a'' not found in workspace.');
end
groups = {'DN_ONSus_RF_UV','DN_ONSus_RF_GRN'};  % edit/extend as needed
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

    % ---------- Fit DoG-CDF to mean curve ----------
    [s_fit, y_fit, p_struct, fit_stats] = fit_DoG_half_CDF(sorted_sizes(:), mu(:));
    results.(gname).fit.s_fit = s_fit;
    results.(gname).fit.y_fit = y_fit;
    results.(gname).fit.params = p_struct;
    results.(gname).fit.stats = fit_stats;

    % Overlay fitted curve
    plot(s_fit, y_fit, '-', 'LineWidth', 2.0, 'Color', colororder(g,:));
    
    % Display fit parameters for mean
    fprintf('\n%s DoG-CDF Fit Results (Mean):\n', gname);
    fprintf('  Shared mean: mu=%.2f\n', p_struct.mu);
    fprintf('  Center sigma=%.2f, Surround sigma=%.2f\n', p_struct.sigma_c, p_struct.sigma_s);
    fprintf('  Model: gain=%.2f, w_s=%.3f, global_bias=%.2f\n', p_struct.gain, p_struct.w_s, p_struct.global_bias);
    fprintf('  Fit quality: loss=%.4f, R²=%.3f\n', fit_stats.fval, fit_stats.r_squared);

    % ---------- Fit DoG-CDF to each individual cell ----------
    fprintf('\nFitting individual cells for group %s...\n', gname);
    
    % Initialize storage for individual cell fits
    results.(gname).individual_fits = struct();
    cell_params = struct();
    param_fields = {'mu', 'sigma_c', 'sigma_s', 'w_s', 'gain', 'global_bias'};
    for pf = 1:length(param_fields)
        cell_params.(param_fields{pf}) = nan(nCells, 1);
    end
    cell_r_squared = nan(nCells, 1);
    cell_rmse = nan(nCells, 1);
    
    for c = 1:nCells
        cell_response = S(:, c);  % Response of cell c to different sizes
        
        % Check if cell has enough valid data points
        valid_pts = ~isnan(cell_response);
        if sum(valid_pts) < 5  % Need at least 5 points for 6 parameters
            fprintf('  Cell %d: Insufficient data points (%d), skipping\n', c, sum(valid_pts));
            continue;
        end
        
        try
            % Fit DoG-CDF to this individual cell
            [s_fit_cell, y_fit_cell, p_struct_cell, fit_stats_cell] = ...
                fit_DoG_half_CDF(sorted_sizes(valid_pts), cell_response(valid_pts));
            
            % Store individual cell fit results
            results.(gname).individual_fits.(sprintf('cell_%d', c)).s_fit = s_fit_cell;
            results.(gname).individual_fits.(sprintf('cell_%d', c)).y_fit = y_fit_cell;
            results.(gname).individual_fits.(sprintf('cell_%d', c)).params = p_struct_cell;
            results.(gname).individual_fits.(sprintf('cell_%d', c)).stats = fit_stats_cell;
            
            % Store parameters for summary statistics
            for pf = 1:length(param_fields)
                cell_params.(param_fields{pf})(c) = p_struct_cell.(param_fields{pf});
            end
            cell_r_squared(c) = fit_stats_cell.r_squared;
            cell_rmse(c) = fit_stats_cell.rmse;
            
            % Create and save individual cell plot
            fig_cell = figure('Color', 'w', 'Position', [100, 100, 600, 400]);
            hold on;
            
            % Plot data points
            plot(sorted_sizes, cell_response, 'o', 'MarkerSize', 8, ...
                'Color', colororder(g,:), 'MarkerFaceColor', colororder(g,:), ...
                'LineWidth', 1.5);
            % Plot fitted curve
            plot(s_fit_cell, y_fit_cell, '-', 'LineWidth', 2.5, 'Color', colororder(g,:));
            % Styling
            grid on;
            xlim([0, max(sorted_sizes)*1.05]);
            xlabel('Spot Diameter');
            ylabel('Mean Firing Rate (Hz)');
            title(sprintf('%s - Cell %d (R²=%.3f)', strrep(gname, '_', '\_'), c, fit_stats_cell.r_squared));
            % Add parameter text
            param_text = sprintf(['Fit Parameters:\n' ...
                'Shared μ=%.2f\n' ...
                'Center σ=%.2f, Surround σ=%.2f\n' ...
                'Weight: w_s=%.3f\n' ...
                'Gain: %.2f, Bias: %.2f'], ...
                p_struct_cell.mu, p_struct_cell.sigma_c, ...
                p_struct_cell.sigma_s, ...
                p_struct_cell.w_s, p_struct_cell.gain, p_struct_cell.global_bias);
            text(0.02, 0.98, param_text, 'Units', 'normalized', ...
                'VerticalAlignment', 'top', 'FontSize', 9, ...
                'BackgroundColor', [1 1 1 0.8], 'EdgeColor', 'k');
            % Save figure
            filename = sprintf('%s_Cell_%d_SpotSizeFit.png', gname, c);
            saveas(fig_cell, fullfile(save_fig_folder, filename));
            close(fig_cell);
            
            fprintf('  Cell %d: R²=%.3f, RMSE=%.3f, saved as %s\n', ...
                c, fit_stats_cell.r_squared, fit_stats_cell.rmse, filename);
            
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
title('Size tuning: mean ± SEM across cells with DoG-CDF fit');
legend(groups, 'Interpreter', 'none', 'Location', 'best');
hold off;

% Save the summary figure
summary_filename = 'SpotSizeAnalysis_Summary.png';
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
            param_fields = {'mu', 'sigma_c', 'sigma_s', 'w_s', 'gain', 'global_bias'};
            param_names = {'Shared μ', 'Center σ', 'Surround σ', 'Surround weight', 'Gain', 'Bias'};
            for pf = 1:length(param_fields)
                if isfield(results.(gname).cell_params_summary, param_fields{pf})
                    param_stats = results.(gname).cell_params_summary.(param_fields{pf});
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