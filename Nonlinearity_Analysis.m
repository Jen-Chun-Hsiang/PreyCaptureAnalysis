%% Derive rec_thr and softness from NL_params (no refit)
% NL_params columns: pfit(:)' center_x center_y
% For 'cdf' fit we assume pfit = [A, mu, sigma, offset]
num_cells = size(NL_params,1);
rec_thr = nan(num_cells,1);
softness = nan(num_cells,1);
for iC = 1:num_cells
    % center_x is NL_params(iC,5)
    rec_thr(iC) = NL_params(iC,5);
    sigma = NL_params(iC,3);
    % Map sigma (cdf width) to softplus softness using width-matching constant
    softness(iC) = max(1e-6, 0.5833 * sigma);
end

% Save these derived params
if ~exist(NL_save_folder, 'dir')
    mkdir(NL_save_folder);
end
save(fullfile(NL_save_folder, 'NL_derived_softplus_params.mat'), 'rec_thr', 'softness', 'NL_params');


fprintf('Average rec_thr: %f\n', mean(rec_thr(~isoutlier(rec_thr))));
%% Grouping: Temporal vs Nasal and ON vs OFF (heuristics)
isTemporal = false(num_cells,1);
isON = false(num_cells,1);
for iC = 1:num_cells
    locstr = '';
    if exist('location','var') && length(location) >= iC && ischar(location{iC})
        locstr = lower(location{iC});
    end
    ctstr = '';
    if exist('cell_type','var') && length(cell_type) >= iC && ischar(cell_type{iC})
        ctstr = lower(cell_type{iC});
    end
    if ~isempty(locstr) && (contains(locstr,'temporal') || contains(locstr,'temp'))
        isTemporal(iC) = true;
    elseif ~isempty(locstr) && (contains(locstr,'nasal') || contains(locstr,'nas'))
        isTemporal(iC) = false;
    else
        if exist('data_sets','var') && length(data_sets) >= iC && ischar(data_sets{iC}) && contains(lower(data_sets{iC}),'temp')
            isTemporal(iC) = true;
        end
    end
    if ~isempty(ctstr) && contains(ctstr,'on')
        isON(iC) = true;
    elseif ~isempty(ctstr) && contains(ctstr,'off')
        isON(iC) = false;
    else
        if exist('cell_type_numeric','var') && length(cell_type_numeric) >= iC
            isON(iC) = cell_type_numeric(iC) == 1;
        elseif exist('data_sets','var') && length(data_sets) >= iC && ischar(data_sets{iC}) && contains(lower(data_sets{iC}),'on')
            isON(iC) = true;
        end
    end
end

%% Stats and comparisons
groups = struct();
for state = 0:1
    for loc = 0:1
        mask = (isON == logical(state)) & (isTemporal == logical(loc));
        groups(state+1,loc+1).mask = mask;
        groups(state+1,loc+1).n = sum(mask);
        groups(state+1,loc+1).rec_thr = rec_thr(mask);
        groups(state+1,loc+1).softness = softness(mask);
    end
end

results = struct();
for state = 0:1
    a_thr = groups(state+1,1).rec_thr; % nasal
    b_thr = groups(state+1,2).rec_thr; % temporal
    a_soft = groups(state+1,1).softness;
    b_soft = groups(state+1,2).softness;
    % t-tests
    try [h1,p1] = ttest2(a_thr,b_thr); catch, h1=NaN; p1=NaN; end
    try p1_rs = ranksum(a_thr,b_thr); catch, p1_rs=NaN; end
    try [h2,p2] = ttest2(a_soft,b_soft); catch, h2=NaN; p2=NaN; end
    try p2_rs = ranksum(a_soft,b_soft); catch, p2_rs=NaN; end
    % effect sizes (Cohen's d)
    d_thr = cohens_d(a_thr,b_thr);
    d_soft = cohens_d(a_soft,b_soft);
    results(state+1).rec_thr = struct('t_p',p1,'ranksum_p',p1_rs,'cohens_d',d_thr);
    results(state+1).softness = struct('t_p',p2,'ranksum_p',p2_rs,'cohens_d',d_soft);
end

%% Compare ON vs OFF combining nasal and temporal
fprintf('\n==== ON vs OFF (combined nasal+temporal) ===\n');
maskON = isON;
maskOFF = ~isON;
on_rec = rec_thr(maskON);
off_rec = rec_thr(maskOFF);
on_soft = softness(maskON);
off_soft = softness(maskOFF);
% tests
try [h_rec,p_rec] = ttest2(on_rec, off_rec); catch, h_rec=NaN; p_rec=NaN; end
try p_rec_rs = ranksum(on_rec, off_rec); catch, p_rec_rs = NaN; end
try [h_soft,p_soft] = ttest2(on_soft, off_soft); catch, h_soft=NaN; p_soft=NaN; end
try p_soft_rs = ranksum(on_soft, off_soft); catch, p_soft_rs = NaN; end
% effect sizes
d_rec = cohens_d(on_rec, off_rec);
d_soft = cohens_d(on_soft, off_soft);
% print
fprintf('\nrec_thr (ON vs OFF):\n');
fprintf('  ON:  mean=%.4g, sd=%.4g, n=%d\n', mean_or_nan(on_rec), std_or_nan(on_rec), sum(maskON));
fprintf('  OFF: mean=%.4g, sd=%.4g, n=%d\n', mean_or_nan(off_rec), std_or_nan(off_rec), sum(maskOFF));
fprintf('  t-test p = %.4g, ranksum p = %.4g, cohens d = %.4g\n', p_rec, p_rec_rs, d_rec);

fprintf('\nsoftness (ON vs OFF):\n');
fprintf('  ON:  mean=%.4g, sd=%.4g, n=%d\n', mean_or_nan(on_soft), std_or_nan(on_soft), sum(maskON));
fprintf('  OFF: mean=%.4g, sd=%.4g, n=%d\n', mean_or_nan(off_soft), std_or_nan(off_soft), sum(maskOFF));
fprintf('  t-test p = %.4g, ranksum p = %.4g, cohens d = %.4g\n', p_soft, p_soft_rs, d_soft);

% Save into results
% Use deal to safely assign ONvsOFF into results struct (handles struct-array cases)
[results.ONvsOFF] = deal(struct('rec_thr', struct('t_p', p_rec, 'ranksum_p', p_rec_rs, 'cohens_d', d_rec), 'softness', struct('t_p', p_soft, 'ranksum_p', p_soft_rs, 'cohens_d', d_soft)));
save(fullfile(NL_save_folder,'NL_derived_stats.mat'),'groups','results','rec_thr','softness');

%% Print summary stats to MATLAB terminal
fprintf('\nNonlinearity derived stats (rec_thr and softness)\n');
for state = 0:1
    state_name = ternary(state, 'ON', 'OFF');
    fprintf('\n---- %s cells ----\n', state_name);
    a_thr = groups(state+1,1).rec_thr; b_thr = groups(state+1,2).rec_thr; % nasal vs temporal
    a_soft = groups(state+1,1).softness; b_soft = groups(state+1,2).softness;
    fprintf('Group counts (Nasal, Temporal): %d, %d\n', groups(state+1,1).n, groups(state+1,2).n);
    fprintf('\nrec_thr:\n');
    fprintf('  Nasal:    mean=%.4g, sd=%.4g\n', mean_or_nan(a_thr), std_or_nan(a_thr));
    fprintf('  Temporal: mean=%.4g, sd=%.4g\n', mean_or_nan(b_thr), std_or_nan(b_thr));
    fprintf('  t-test p = %.4g, ranksum p = %.4g, cohens d = %.4g\n', results(state+1).rec_thr.t_p, results(state+1).rec_thr.ranksum_p, results(state+1).rec_thr.cohens_d);

    fprintf('\nsoftness:\n');
    fprintf('  Nasal:    mean=%.4g, sd=%.4g\n', mean_or_nan(a_soft), std_or_nan(a_soft));
    fprintf('  Temporal: mean=%.4g, sd=%.4g\n', mean_or_nan(b_soft), std_or_nan(b_soft));
    fprintf('  t-test p = %.4g, ranksum p = %.4g, cohens d = %.4g\n', results(state+1).softness.t_p, results(state+1).softness.ranksum_p, results(state+1).softness.cohens_d);
end

%% Plot boxplots
figure; clf; tiledlayout(1,2);
nexttile; hold on;
boxdata = [];
glabels = {};
for state=0:1
    for loc=0:1
        vals = groups(state+1,loc+1).rec_thr;
        boxdata = [boxdata; vals(:)];
        glabels = [glabels; repmat({sprintf('%s-%s', ternary(state,'OFF','ON'), ternary(loc,'Nasal','Temporal'))}, length(vals),1)];
    end
end
boxplot(boxdata, glabels, 'LabelOrientation','inline'); title('rec\_thr'); ylabel('rec\_thr');
nexttile; hold on;
boxdata = [];
glabels = {};
for state=0:1
    for loc=0:1
        vals = groups(state+1,loc+1).softness;
        boxdata = [boxdata; vals(:)];
        glabels = [glabels; repmat({sprintf('%s-%s', ternary(state,'OFF','ON'), ternary(loc,'Nasal','Temporal'))}, length(vals),1)];
    end
end
boxplot(boxdata, glabels, 'LabelOrientation','inline'); title('softness'); ylabel('softness');
saveas(gcf, fullfile(NL_save_folder,'NL_derived_boxplots.png'));

%% Normality-checked comparisons: use t-test only if both groups pass Lilliefors normality test, else use ranksum
fprintf('\n==== Normality-checked comparisons (Lilliefors) ===\n');
alpha = 0.05;
results.normality_checked = struct();
for state = 0:1
    a_thr = groups(state+1,1).rec_thr; b_thr = groups(state+1,2).rec_thr;
    a_soft = groups(state+1,1).softness; b_soft = groups(state+1,2).softness;
    res_thr = compare_groups(a_thr,b_thr,alpha);
    res_soft = compare_groups(a_soft,b_soft,alpha);
    results.normality_checked(state+1).rec_thr = res_thr;
    results.normality_checked(state+1).softness = res_soft;
    % print
    state_name = ternary(state,'ON','OFF');
    fprintf('\n-- %s (Nasal vs Temporal) --\n', state_name);
    fprintf('rec_thr: method=%s, p=%.4g, ranksum_p=%.4g, cohens_d=%.4g\n', res_thr.method, res_thr.p, res_thr.ranksum_p, res_thr.cohens_d);
    fprintf('softness: method=%s, p=%.4g, ranksum_p=%.4g, cohens_d=%.4g\n', res_soft.method, res_soft.p, res_soft.ranksum_p, res_soft.cohens_d);
end

% ON vs OFF combined (normality-checked)
res_onoff_thr = compare_groups(rec_thr(isON), rec_thr(~isON), alpha);
res_onoff_soft = compare_groups(softness(isON), softness(~isON), alpha);
results.normality_checked.ONvsOFF = struct('rec_thr', res_onoff_thr, 'softness', res_onoff_soft);
fprintf('\n-- ON vs OFF (combined) --\n');
fprintf('rec_thr: method=%s, p=%.4g, ranksum_p=%.4g, cohens_d=%.4g\n', res_onoff_thr.method, res_onoff_thr.p, res_onoff_thr.ranksum_p, res_onoff_thr.cohens_d);
fprintf('softness: method=%s, p=%.4g, ranksum_p=%.4g, cohens_d=%.4g\n', res_onoff_soft.method, res_onoff_soft.p, res_onoff_soft.ranksum_p, res_onoff_soft.cohens_d);

% Save updated results
save(fullfile(NL_save_folder,'NL_derived_stats.mat'),'groups','results','rec_thr','softness');

%% Helper functions
function d = cohens_d(x,y)
    x = x(:); y = y(:);
    nx = numel(x); ny = numel(y);
    if isempty(x) || isempty(y)
        d = NaN; return
    end
    sx = var(x,1); sy = var(y,1);
    s_pooled = sqrt(((nx-1)*sx + (ny-1)*sy) / (nx+ny-2));
    if s_pooled == 0
        d = NaN; return
    end
    d = (mean(x) - mean(y)) / s_pooled;
end

function out = ternary(cond,a,b)
    if cond
        out = a;
    else
        out = b;
    end
end

%% small helpers for printing
function m = mean_or_nan(x)
    if isempty(x)
        m = NaN;
    else
        m = mean(x);
    end
end

function s = std_or_nan(x)
    if isempty(x)
        s = NaN;
    else
        s = std(x);
    end
end

%% Helper: compare two groups with normality check
function out = compare_groups(x,y,alpha)
    if nargin < 3, alpha = 0.05; end
    x = x(:); y = y(:);
    x = x(~isnan(x)); y = y(~isnan(y));
    out.nx = numel(x); out.ny = numel(y);
    out.meanx = mean_or_nan(x); out.meany = mean_or_nan(y);
    out.sdx = std_or_nan(x); out.sdy = std_or_nan(y);
    % check normality with Lilliefors (lillietest). Require at least 3 samples each.
    if out.nx >= 3 && out.ny >= 3
        try
            hx = lillietest(x,'Alpha',alpha);
            hy = lillietest(y,'Alpha',alpha);
        catch
            hx = 1; hy = 1; % if test fails, treat as non-normal
        end
    else
        hx = 1; hy = 1; % small sample -> non-normal
    end
    out.normal = (~hx) && (~hy);
    if out.normal
        % use Welch t-test
        try
            [~,pt] = ttest2(x,y,'Vartype','unequal');
        catch
            pt = NaN;
        end
        out.method = 'ttest2';
        out.p = pt;
    else
        try
            pr = ranksum(x,y);
        catch
            pr = NaN;
        end
        out.method = 'ranksum';
        out.p = pr;
    end
    % always compute both for reference
    try
        pr_all = ranksum(x,y);
    catch
        pr_all = NaN;
    end
    try
        [~,pt_all] = ttest2(x,y,'Vartype','unequal');
    catch
        pt_all = NaN;
    end
    out.ranksum_p = pr_all;
    out.ttest_p = pt_all;
    out.cohens_d = cohens_d(x,y);
end
