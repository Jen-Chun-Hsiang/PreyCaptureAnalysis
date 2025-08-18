
group_id = 2;
switch group_id
    case 1
        % ONSustained AcuteZone vs ONSustained DN
        x = ONSus_AcuteZone_SI;
        y = ONSus_DN_SI;
        group_names = {'ONSus AcuteZone', 'ONSus DN'};
    case 2
        % OFFSustained AcuteZone vs OFFSustained DN
        x = OFFSus_AcuteZone_SI;
        y = OFFSus_DN_SI;
        group_names = {'OFFSus AcuteZone', 'OFFSus DN'};
    otherwise
        error('Invalid group_id');
end

[hx, px] = lillietest(x);
[hy, py] = lillietest(y);

if hx == 0 && hy == 0
    disp('Both groups are likely normal. ttest2 is appropriate.');
    [h, p] = ttest2(x, y);
    fprintf('ttest2: p = %.4g, h = %d\n', p, h);
else
    disp('At least one group is not normal. Consider using ranksum instead.');
    [h, p] = kstest2(x, y);
    fprintf('kstest2: p = %.4g, h = %d\n', p, h);
end

% Test if x is significantly different from zero
if hx == 0
    [hx0, px0] = ttest(x);
    fprintf('x: t-test vs zero, p = %.4g\n', px0);
else
    [px0, hx0] = signrank(x);
    fprintf('x: signrank vs zero, p = %.4g\n', px0);
end

% Test if y is significantly different from zero
if hy == 0
    [hy0, py0] = ttest(y);
    fprintf('y: t-test vs zero, p = %.4g\n', py0);
else
    [py0, hy0] = signrank(y);
    fprintf('y: signrank vs zero, p = %.4g\n', py0);
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
ylabel('Value');
title('Group Comparison with Individual Data Points');
box off;
hold
