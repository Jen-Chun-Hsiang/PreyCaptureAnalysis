%% 
% Run  WhiteNoise_ONOFFalpha_Comparison.m first
% Improved bar plot with matching-colored error bars, white edges, and no caps
figure;
hold on

% -------------------------------
% 1. Calculate normalized peaks
% -------------------------------

peak_nasal    = biphase(locs==0);
peak_temporal = biphase(locs==1);

% Compute group means and SEMs
means = [ mean(peak_nasal),  mean(peak_temporal) ];
sems  = [ std(peak_nasal)/sqrt(numel(peak_nasal)),  std(peak_temporal)/sqrt(numel(peak_temporal)) ];

% -------------------------------
% 2. Define colors
% -------------------------------
% Choose a distinct color for each bar (RGB triples in [0,1])
colors = [
    26 51 233;   % Nasal – a medium blue
    221 52 137;    % Temporal – a brick red
]/255;

% -------------------------------
% 3. Create the bar plot
% -------------------------------
% Use 'FaceColor','flat' so we can assign each bar its own color via CData
b = bar(1:2, means, ...
    'FaceColor', 'flat', ...
    'EdgeColor', 'white', ...
    'LineWidth', 1.5);  % White edges, slightly thicker

% Assign each bar its corresponding RGB color
for ii = 1:2
    b.CData(ii, :) = colors(ii, :);
end

% -------------------------------
% 4. Add error bars (no caps)
% -------------------------------
% Plot each error bar separately to match the bar color exactly
for ii = 1:2
    errorbar( ...
        ii, ...                 % x-position
        means(ii), ...          % y-value (mean)
        sems(ii), ...           % negative error
        sems(ii), ...           % positive error
        'LineStyle',   'none', ...
        'LineWidth',    1.5, ...
        'Color',        0.2*ones(1, 3), ...
        'CapSize',      0 ...   % No caps on the error bars
    );
end

% -------------------------------
% 5. Final Aesthetics
% -------------------------------
xlim([0.5, 2.5]);                % Tighten x-limits around the bars
xticks([1 2]);
xticklabels({sprintf('Nasal (n=%d)', numel(peak_nasal)),...
    sprintf('Temporal (n=%d)', numel(peak_temporal))});
ylabel('Biphasic','FontSize', 12);
ylim([0 0.35]);
yticks(0:0.1:0.3);
yticklabels({'0', '0.1', '0.2', '0.3'})
set(gca, ...
    'Box',      'off',           ...  % Remove top and right box lines
    'FontSize', 12,             ...  % Axis font size
    'LineWidth', 1);                 % Axis line width

% Optional: Add gridlines if desired
% grid on
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

hold off
