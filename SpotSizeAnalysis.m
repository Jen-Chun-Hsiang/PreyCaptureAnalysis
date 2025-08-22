% Analyze and plot size tuning (mean during stimulus ± SEM across cells)
% Requires: struct 'a' in workspace with fields:
%   - x1 (T x 1), y_labels (1 x S), and group arrays (T x (S * Ncells))
% Example groups: 'DN_ONSus_RF_UV', 'DN_ONSus_RF_GRN'

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
end

% ---------------- figure styling ----------------
grid on;
xlabel('Diameter');
ylabel('Mean firing rate during stimulus (Hz)');
title('Size tuning: mean ± SEM across cells');
legend(groups, 'Interpreter', 'none', 'Location', 'best');
hold off;

% Results struct is left in workspace as "results"