
group_id = 8;
clear measured_y
switch group_id
    case 1
        % ONSustained AcuteZone vs ONSustained DN
        x = ONSus_AcuteZone(:, 1);
        y = ONSus_DorsalNasal(:, 1);
        measured_y = '(S-C)/C';
        group_names = {'ONSus AcuteZone', 'ONSus DN'};
    case 2
        % OFFSustained AcuteZone vs OFFSustained DN
        x = OFFSus_AcuteZone(:, 1);
        y = OFFSus_DorsalNasal(:, 1);
        measured_y = '(S-C)/C';
        group_names = {'OFFSus AcuteZone', 'OFFSus DN'};
    case 3
        x = -TF_time2peak(Ids{3});
        y = -TF_time2peak(Ids{4});
        measured_y = 'Time to peak (ms)';
        group_names = {'ONSus AcuteZone', 'ONSus DN'};
    case 4
        x = -TF_time2peak(Ids{5});
        y = -TF_time2peak(Ids{6});
        measured_y = 'Time to peak (ms)';
        group_names = {'OFFSus AcuteZone', 'OFFSus DN'};
    case 5
        rf_area = rf_pixels*4.375^2;
        x = rf_area(Ids{3});
        y = rf_area(Ids{4});
        measured_y = 'RF Area (um^2)';
        group_names = {'ONSus AcuteZone', 'ONSus DN'};
    case 6
        rf_area = rf_pixels*4.375^2;
        x = rf_area(Ids{5});
        y = rf_area(Ids{6});
        measured_y = 'RF Area (um^2)';
        group_names = {'OFFSus AcuteZone', 'OFFSus DN'};
    case 7
        x = avg_rad(Ids{3});
        y = avg_rad(Ids{4});
        measured_y = 'Average diameter (um)';
        group_names = {'ONSus AcuteZone', 'ONSus DN'};
    case 8
        x = avg_rad(Ids{5});
        y = avg_rad(Ids{6});
        measured_y = 'Average diameter (um)';
        group_names = {'OFFSus AcuteZone', 'OFFSus DN'};
    case 9
        x = TF_biphasic_peaks(Ids{3});
        y = TF_biphasic_peaks(Ids{4});
        measured_y = 'Biphasic';
        group_names = {'ONSus AcuteZone', 'ONSus DN'};
    case 10
        x = TF_biphasic_peaks(Ids{5});
        y = TF_biphasic_peaks(Ids{6});
        measured_y = 'Biphasic';
        group_names = {'OFFSus AcuteZone', 'OFFSus DN'};
    otherwise
        error('Invalid group_id');
end
x = x(~isnan(x));
y = y(~isnan(y));

%%
[hx, px] = lillietest(x);
[hy, py] = lillietest(y);

if hx == 0 && hy == 0
    disp('Both groups are likely normal. ttest2 is appropriate.');
    [h, p, ci, stats] = ttest2(x, y);
    fprintf('ttest2: p = %.4g, h = %d\n', p, h);
    
    % Calculate Cohen's d effect size
    pooled_std = sqrt(((length(x)-1)*var(x) + (length(y)-1)*var(y)) / (length(x) + length(y) - 2));
    cohens_d = (mean(x) - mean(y)) / pooled_std;
    fprintf('Cohen''s d = %.3f\n', cohens_d);
    
    % Interpret effect size
    if abs(cohens_d) < 0.2
        effect_interp = 'negligible';
    elseif abs(cohens_d) < 0.5
        effect_interp = 'small';
    elseif abs(cohens_d) < 0.8
        effect_interp = 'medium';
    else
        effect_interp = 'large';
    end
    fprintf('Effect size interpretation: %s\n', effect_interp);
    
else
    disp('At least one group is not normal. Consider using ranksum instead.');
    [p, h, stats] = ranksum(x, y);
    fprintf('ranksum: p = %.4g, h = %d\n', p, h);
    
    % Calculate effect size r for Mann-Whitney U test
    % For ranksum, calculate z-score from ranksum statistic
    n1 = length(x);
    n2 = length(y);
    U = stats.ranksum - n1*(n1+1)/2;  % Convert ranksum to U statistic
    mu_U = n1*n2/2;
    sigma_U = sqrt(n1*n2*(n1+n2+1)/12);
    z_score = (U - mu_U) / sigma_U;
    
    n_total = n1 + n2;
    effect_r = abs(z_score) / sqrt(n_total);
    fprintf('Effect size r = %.3f\n', effect_r);
    
    % Interpret effect size
    if effect_r < 0.1
        effect_interp = 'negligible';
    elseif effect_r < 0.3
        effect_interp = 'small';
    elseif effect_r < 0.5
        effect_interp = 'medium';
    else
        effect_interp = 'large';
    end
    fprintf('Effect size interpretation: %s\n', effect_interp);
end

% Test if x is significantly different from zero
if hx == 0
    [hx0, px0, ci, stats] = ttest(x);
    fprintf('x: t-test vs zero, p = %.4g\n', px0);
    % Cohen's d for one-sample t-test
    cohens_d_x = mean(x) / std(x);
    fprintf('x: Cohen''s d vs zero = %.3f\n', cohens_d_x);
else
    [px0, hx0, stats] = signrank(x);
    fprintf('x: signrank vs zero, p = %.4g\n', px0);
    % Effect size r for signed rank test
    % Calculate z-score manually for signrank
    n = length(x);
    W = stats.signedrank;
    mu_W = n*(n+1)/4;
    sigma_W = sqrt(n*(n+1)*(2*n+1)/24);
    z_score_x = (W - mu_W) / sigma_W;
    effect_r_x = abs(z_score_x) / sqrt(n);
    fprintf('x: Effect size r vs zero = %.3f\n', effect_r_x);
end

% Test if y is significantly different from zero
if hy == 0
    [hy0, py0, ci, stats] = ttest(y);
    fprintf('y: t-test vs zero, p = %.4g\n', py0);
    % Cohen's d for one-sample t-test
    cohens_d_y = mean(y) / std(y);
    fprintf('y: Cohen''s d vs zero = %.3f\n', cohens_d_y);
else
    [py0, hy0, stats] = signrank(y);
    fprintf('y: signrank vs zero, p = %.4g\n', py0);
    % Effect size r for signed rank test
    % Calculate z-score manually for signrank
    n = length(y);
    W = stats.signedrank;
    mu_W = n*(n+1)/4;
    sigma_W = sqrt(n*(n+1)*(2*n+1)/24);
    z_score_y = (W - mu_W) / sigma_W;
    effect_r_y = abs(z_score_y) / sqrt(n);
    fprintf('y: Effect size r vs zero = %.3f\n', effect_r_y);
end

% filepath: \\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\rf_property_stats.m
% ...existing code...
% After statistical tests

figure;
hold on;

% Data for plotting
data = {x, y};

% Bar chart of means
bar_vals = cellfun(@mean, data);
bar(1:2, bar_vals, 'FaceColor', 'flat', 'CData', [0.2 0.6 1; 1 0.4 0.4]);

% Overlay individual data points (jittered)
for i = 1:2
    jitter = (rand(size(data{i})) - 0.5) * 0.15;
    scatter(i + jitter, data{i}, 40, 'k', 'filled', 'MarkerFaceAlpha', 0.7);
end

set(gca, 'XTick', 1:2, 'XTickLabel', group_names);
ylabel(sprintf('%s', measured_y));
title(sprintf('Comparison of %s (n=%d) and %s (n=%d) p-val=%.4g \n', group_names{1}, length(x), group_names{2}, length(y), p));
box off;
hold off;
