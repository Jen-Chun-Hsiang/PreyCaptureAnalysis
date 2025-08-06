clear; close all; clc;
recordingnames = {'f081224', 'c081224', 'd081024', 'c081024', 'a081024'};
loc_ids = [1 1 2 2 2];
num_rec = length(recordingnames);

load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\StationarySpot\';

%%

for i = 1:num_rec
    recordingname = recordingnames{i};
    saveFileName = sprintf('%s_stationary_spot.mat', recordingname);
    load(sprintf('%s/%s', load_data_folder, saveFileName), 'diameters', 'ct', 'contrast_name','Data');
    
    if i ==1
        peaks = nan(num_rec, numel(diameters));
    end
    tids = find(ct > 0.58 & ct < 0.8);
    peaks(i, :) = max(mean(Data(:, 2, :, tids), 3), [], 4);
end

colors = [1, 0, 0; 0, 1, 0];
figure; hold on
for i = 1:size(peaks, 1)
    plot(diameters, peaks(i, :), 'Color', colors(loc_ids(i), :));
end

figure; hold on
peak_nasal = mean(peaks(loc_ids==1, end), 2)./mean(peaks(loc_ids==1, 4:5), 2);
peak_temporal = mean(peaks(loc_ids==2, end), 2)./mean(peaks(loc_ids==2, 4:5), 2);
bar([1, 2], [mean(peak_nasal) mean(peak_temporal)])
errorbar([1, 2], [mean(peak_nasal) mean(peak_temporal)],...
    [std(peak_nasal)/sqrt(2), std(peak_temporal)/sqrt(3)]);

%%
% Improved bar plot with matching-colored error bars, white edges, and no caps
figure;
hold on

% -------------------------------
% 1. Calculate normalized peaks
% -------------------------------
% For 'nasal' (loc_ids == 1):
nasal_vals    = peaks(loc_ids == 1, :);  % Extract rows where loc_ids == 1
peak_nasal    = mean(nasal_vals(:, end), 2) ./ mean(nasal_vals(:, 4:5), 2);

% For 'temporal' (loc_ids == 2):
temporal_vals = peaks(loc_ids == 2, :);  % Extract rows where loc_ids == 2
peak_temporal = mean(temporal_vals(:, end), 2) ./ mean(temporal_vals(:, 4:5), 2);

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
xticklabels({'Nasal (n=2)','Temporal (n=3)'});
ylabel('Large vs peak response ratio','FontSize', 12);
ylim([0.6 1.2]);
yticks(0.6:0.2:1.2);
yticklabels({'0.6', '0.8', '1.0', '1.2'})
set(gca, ...
    'Box',      'off',           ...  % Remove top and right box lines
    'FontSize', 12,             ...  % Axis font size
    'LineWidth', 1);                 % Axis line width

% Optional: Add gridlines if desired
% grid on
% set(gca, 'GridLineStyle', ':', 'GridAlpha', 0.3);

hold off
