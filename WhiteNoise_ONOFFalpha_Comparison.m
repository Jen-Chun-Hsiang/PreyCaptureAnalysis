close all; clear; clc;
%%
data_sets = {'e100724', 'f100724', 'a101224', 'b101224',  'c101224',  'd101224', 'e101224',...
             'b101424', 'c101424', 'd101424', 'e101424',  'a101624',  'b101624', 'd101624', 'e101624',...
             'b101924', 'c101924', 'd101924', 'e101924',  'b103124',  'e103124', 'a110424', 'b110424',...
             'c110424', 'd110424', 'e110424', 'f110424',  'g110424',  'a110924', 'b110924', 'c110924',...
             'a111224'};
cell_type = {'OFF',      'OFF',    'OFF',      'ON',       'OFF',     'ON',      'OFF',...
             'OFF',      'OFF',    'ON',       'ON',       'ON',      'ON',      'ON',       'ON',...
             'ON',       'OFF',    'OFF',      'OFF',      'ON',      'OFF',     'ON',       'ON',...
             'ON',       'ON',     'OFF',      'ON',       'OFF',     'ON',      'OFF',      'OFF',...
             'ON'};
location =  {'Temporal', 'Temporal','Nasal',   'Nasal',    'Nasal',   'Nasal',   'Nasal',...
             'Temporal', 'Temporal','Temporal','Temporal', 'Nasal',   'Nasal',   'Nasal',    'Nasal',...
             'Nasal',    'Nasal',   'Nasal',   'Nasal',    'Temporal','Temporal','Temporal', 'Temporal',...
             'Temporal', 'Temporal','Temporal','Temporal', 'Temporal','Temporal','Temporal', 'Temporal',...
             'Temporal',...
             };

process_version = 'GaussianFitting_processed_081425_1.mat';
is_normalized_tf = 1;
save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
%% Check files
% Load Excel sheet
excelFile = 'PatchRecordingAlphaRGCs_PreyCapture.xlsx'; % <-- change to your actual file name
T = readtable(excelFile, 'Sheet', 'Sheet1');

% Convert Excel date to string for matching
%excelDayStr = cellfun(@(x) datestr(x, 'mmddyy'), T.Day, 'UniformOutput', false);
excelDayStr = arrayfun(@(x) datestr(x, 'mmddyy'), T.Day, 'UniformOutput', false);

for i = 1:length(data_sets)
    fprintf('Checking %s...%d \n', data_sets{i}, i);
    ds = data_sets{i};
    % Extract letter and date from data_sets
    cellLetter = ds(1);
    dateStr = ds(2:end);
    
    % Find matching row in Excel
    matchIdx = find(strcmp(excelDayStr, dateStr) & strcmp(T.Cell, cellLetter));
    if isempty(matchIdx)
        fprintf('No match for %s in Excel sheet!\n', ds);
        error('Stopping due to missing match.');
    end
    
    % Check cell_type
    excelType = T.('ON_OFF'){matchIdx};
    if ~strcmpi(cell_type{i}, excelType)
        error(sprintf('Type mismatch for %s: data_sets=%s, Excel=%s\n', ds, cell_type{i}, excelType));
    end
    
    % Check location
    excelLoc = T.TN_axis{matchIdx};
    expectedLoc = location{i};
    if (strcmpi(expectedLoc, 'Temporal') && ~strcmpi(excelLoc, 'T')) || ...
       (strcmpi(expectedLoc, 'Nasal') && ~strcmpi(excelLoc, 'N'))
        error(printf('Location mismatch for %s: data_sets=%s, Excel=%s\n', ds, expectedLoc, excelLoc));
    end
    fprintf('Checked %s: Type=%s, Location=%s PASS \n', ds, excelType, excelLoc);
end
%%
clear Data
num_set = length(data_sets);
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
for i = 1:num_set
    file_name = sprintf('%s.mat', data_sets{i});
    Data{i} = load(fullfile(folder_name, file_name), 'stdSTA', 'tRF');
end

%%
Fz = 100;
WinT = [-0.5 0];
t = WinT(1):1/Fz:WinT(end);
%% digitalize parameters for cell comparison
numeric_parts = regexp(data_sets, '\d+', 'match');
numeric_values = cellfun(@str2double, numeric_parts);
[date_value, ~, date_ids] = unique(numeric_values);
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), location);
cell_type_numeric = cellfun(@(x) strcmp(x, 'ON'), cell_type);
%%
GetWhiteNoise_Nonlinearity
visualize_nonlinear_curve
keyboard;


%% Skip process to section [AA] for loading data
%% show difference in temporal filter



% Check if processed file exists before running spatial receptive field fitting
processedFile = fullfile(folder_name, process_version);
if exist(processedFile, 'file')
    disp('Processed file found. Loading instead of rerunning fitting.');
    load(processedFile);
else
    ct = t(2:end);
    Trace = nan(num_set,  length(ct));

    save_tf_folder = fullfile(save_folder, 'TF_Fits');
    if ~exist(save_tf_folder, 'dir')
        mkdir(save_tf_folder);
    end
    x = 1:length(ct);
    figure; hold on
    for k = 1:num_set
        clc
        fprintf('Processing %s... %d/%d \n', data_sets{k}, k, num_set);
        csig = Data{k}.tRF;
        if is_normalized_tf
            csig = csig./max(abs(csig));
        end
        OptW = GaussianTemporalFilter2(csig');

        tf = gaussmf(x, [OptW(1) OptW(3)])*OptW(5) - gaussmf(x, [OptW(2) OptW(4)])*(OptW(6)*OptW(5)) + OptW(7);

        if k == 1
            Gauss_TF_est = nan(num_set, length(OptW));
        end
        Gauss_TF_est(k, :) = OptW;
        % switch lower(cell_type{k})
        %     case 'on'
        %         plot(ct, csig, 'Color', [247, 224, 12]/255);
        %     case 'off'
        %         plot(ct, csig, 'Color', 0.4*ones(1, 3));
        % end
        Trace(k, :) = csig;

        figure('Visible','off'); hold on
        plot(ct, csig, 'k', 'LineWidth', 1.5); % normalized csig
        plot(ct, tf, 'r--', 'LineWidth', 1.5); % fitted tf
        legend('Normalized csig', 'Fitted TF');
        xlabel('Time (s)');
        ylabel('Normalized Value');
        title(sprintf('TF Fit: %s (%s, %s)', data_sets{k}, cell_type{k}, location{k}));
        set(gca, 'Box', 'off');
        saveas(gcf, fullfile(save_tf_folder, sprintf('TFfit_%s.png', data_sets{k})));
        close(gcf);
    end
    h1 = plot(ct, squeeze(mean(Trace(cell_type_numeric==1, :), 1)), 'Color', [245 182 66]/255, 'LineWidth', 2);
    h2 = plot(ct, squeeze(mean(Trace(cell_type_numeric==0, :), 1)), 'Color', 0*ones(1, 3), 'LineWidth', 2);
    xlabel('Time (s)');
    xticks(-0.5:0.25:0)
    xticklabels({'-0.5', '0.25', '0'});
    if is_normalized_tf
        yticks(-1:0.5:1)
        yticklabels({'-1', '', '0', '', '1'});
    else
        yticks(-100:50:100)
        yticklabels({'-100', '', '0', '', '100'});
    end
    ylabel('Average stimulus value');
    legend([h1, h2], 'ON', 'OFF');
    keyboard;
end
%%


%%
isON = cell_type_numeric == 1;
% isON = Gauss_TF_est(:, 3) >= Gauss_TF_est(:, 4);
% isON_check = Gauss_TF_est(:, 5) > Gauss_TF_est(:, 6); % more accurate
% assert(sum(isON-isON_check)== 0);
%%
TF_time2peak = nan(num_set, 1);
TF_width = nan(num_set, 1);
TF_biphasic_peaks = nan(num_set, 1);
TF_biphasic_stregth = nan(num_set, 1);
hwith_thr = 0.5;
nt = size(Trace, 2);
for i = 1:num_set
    csig = Data{i}.tRF(:)';
    % Interpolate for finer resolution
    csig_interp = interp1(linspace(0, 1, nt), csig, linspace(0, 1, 1000), 'cubic');
    t_interp = linspace(WinT(1), WinT(end), 1000);
    if isON(i)
        [~, peak_idx] = max(csig_interp);
    else
        [~, peak_idx] = min(csig_interp);
    end
    TF_time2peak(i) = (t_interp(peak_idx))*1000; % ms
    maxv = max(csig_interp);
    minv = abs(min(csig_interp));
    if isON(i)
        csig_thr = csig_interp > hwith_thr;
        TF_biphasic_peaks(i) = 1-2*abs(maxv/(maxv+minv)-0.5);
        TF_biphasic_stregth(i) = Gauss_TF_est(i, 6)./Gauss_TF_est(i, 5);
    else
        csig_thr = csig_interp < -hwith_thr;
        TF_biphasic_peaks(i) = 1-2*abs(minv/(maxv+minv)-0.5);
        TF_biphasic_stregth(i) = Gauss_TF_est(i, 5)./Gauss_TF_est(i, 6);
    end
    TF_width(i) = sum(csig_thr)*(1/Fz)*nt;
end

%% 
Disp_Type = 'ON';
cell_type_id = strcmpi(Disp_Type, 'ON');
type_ids = find(cell_type_numeric == cell_type_id)';
biphase = Gauss_TF_est(type_ids, 5:6);
if cell_type_id == 1
    biphase = biphase(:, 2)./biphase(:, 1);
else
    biphase = biphase(:, 1)./biphase(:, 2);
end
colors = parula(256);
color_ids = round(((biphase-min(biphase))*255/range(biphase))+1);
figure; 
subplot(1, 2, 1); hold on
for k = 1:numel(type_ids)
    csig = Data{type_ids(k)}.tRF;
    if is_normalized_tf
        csig = csig./max(abs(csig));
    end
    plot(ct, csig, 'Color', colors(color_ids(k), :));
end

subplot(1, 2, 2); hold on
for k = 1:num_set
    if strcmpi(cell_type{k}, Disp_Type)
        csig = Data{k}.tRF;
        if is_normalized_tf
            csig = csig./max(abs(csig));
        end
        switch lower(location{k})
            case 'temporal'
                plot(ct, csig, 'Color', [66 182 245]/255);
            case 'nasal'
                plot(ct, csig, 'Color', [247 153 205]/255);
        end
    end
end
h1 = plot(ct, squeeze(mean(Trace(cell_type_numeric==cell_type_id & location_type_numeric == 1, :), 1)),...
    'Color', [27 59 242]/255, 'LineWidth', 2);
h2 = plot(ct, squeeze(mean(Trace(cell_type_numeric==cell_type_id & location_type_numeric == 0, :), 1)),...
    'Color', [242 27 145]/255, 'LineWidth', 2);
xlabel('Time (s)');
xticks(-0.5:0.25:0)
xticklabels({'-0.5', '0.25', '0'});
if is_normalized_tf
    yticks(-1:0.5:1)
    yticklabels({'-1', '', '0', '', '1'});
else
    yticks(-100:50:100)
    yticklabels({'-100', '', '0', '', '100'});
end
ylabel('Average stimulus value');
legend([h1, h2], 'Temporal', 'Nasal')
%%
locs = location_type_numeric(type_ids);
x = biphase(locs==0);
y = biphase(locs==1);
% p = ranksum(biphase(locs==0),biphase(locs==1));
[~, p] = ttest2(x, y)
%%
keyboard

%% [AA] save data given a time stamp of processing

process_version = 'GaussianFitting_processed_082025_1.mat';

% Check if processed file exists before running spatial receptive field fitting
processedFile = fullfile(folder_name, process_version);
if exist(processedFile, 'file')
    disp('Processed file found. Loading instead of rerunning fitting.');
    load(processedFile);
else
    % Run spatial receptive field fitting with bounds
    num_gauss = 1;
    image = Data{1}.stdSTA';
    initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];
    initial_params = repmat(initial_params, num_gauss, 1);
    initial_params(:, 1:2) = initial_params(:, 1:2) + 10*rand(num_gauss, 2);
    initial_params = initial_params(:);

    % Example bounds (adjust as needed)
    lb = [1, 1, 5, 5, -pi, -10, 0.5, 5, 5, 0];   % lower bounds
    ub = [size(image,1), size(image,2), 200, 200, pi, 10, 10, 1000, 1000, 10]; % upper bounds
    lb = repmat(lb, num_gauss, 1);
    ub = repmat(ub, num_gauss, 1);

    options = optimoptions('fmincon', 'Display', 'off', 'MaxFunctionEvaluations', 600*length(initial_params));
    for k = 1:num_set
        image = Data{k}.stdSTA';
        objective_function = @(params) 1-corr(image(:), reshape(gaussian_multi(params, image, num_gauss), [], 1));
        optimal_params = fmincon(objective_function, initial_params, [], [], [], [], lb, ub, [], options);
        if k == 1
            gauss_est = nan(num_set, length(optimal_params));
        end
        gauss_est(k, :) = optimal_params;
        clc
        fprintf('Gaussian fitting... (%d/%d)', k, num_set);
    end
    % Check alignment before saving
    for k = 1:num_set
        assert(~isempty(Data{k}), 'Missing Data at index %d', k);
        assert(size(gauss_est,1) >= k, 'gauss_est missing entry at index %d', k);
        assert(length(cell_type_numeric) >= k, 'cell_type_numeric missing entry at index %d', k);
        assert(length(location_type_numeric) >= k, 'location_type_numeric missing entry at index %d', k);
    end
    save(processedFile, 'gauss_est', 'Gauss_TF_est', 'ct', 'Trace', 'location_type_numeric', 'data_sets',...
        'cell_type_numeric', 'cell_type', 'location');
    % Save split group variables to the same file, appending
    split_gauss_groups
    if exist('gauss_est_ON_temporal', 'var')
        save(processedFile, 'gauss_est_ON_temporal', 'gauss_est_ON_nasal', 'gauss_est_OFF_temporal', 'gauss_est_OFF_nasal', ...
            'Gauss_TF_est_ON_temporal', 'Gauss_TF_est_ON_nasal', 'Gauss_TF_est_OFF_temporal', 'Gauss_TF_est_OFF_nasal', '-append');
    end
end
if exist('NL_curves', 'var')
    save(processedFile, 'NL_curves', 'NL_params', '-append');
end

split_est_parameter_by_celltype


%%

elipse_ratio = min(gauss_est(:, 3:4), [], 2)./max(gauss_est(:, 3:4), [], 2);
% avg_rad = sqrt(gauss_est(:, 3).^2 + gauss_est(:, 4).^2)*1.5*4.375;
avg_rad = 2*sqrt(gauss_est(:, 3).*gauss_est(:, 4))*2*4.375;
surround_center = abs(gauss_est(:, 10)./gauss_est(:, 7));
%%
keyboard;
%% Get area size 
FindThreshold_MeanMinusKStd
%%
clc
Colors = lines(6); %
clear Ids ylims ylab values
Ids{1} = cell_type_numeric == 1;
Ids{2} = cell_type_numeric == 0;
Ids{3} = cell_type_numeric == 1 & location_type_numeric == 1;
Ids{4} = cell_type_numeric == 1 & location_type_numeric == 0;
Ids{5} = cell_type_numeric == 0 & location_type_numeric == 1;
Ids{6} = cell_type_numeric == 0 & location_type_numeric == 0;
barlabels = {'ON', 'OFF', 'ON-temporal', 'ON-nasal', 'OFF-temporal', 'OFF-nasal'};
num_n = cellfun(@sum, Ids);
eval_target = 'area';
switch lower(eval_target)
    case 'nl_baseline'
        zero_id = find(x_sample_range==0);
        values = NL_curves(:, zero_id);
        ylab = 'NL baseline param';
    case 'nl_baseline_param'
        values = NL_params(:, 4);
        ylims = [-0.2 25];
        ytick = 0:10:20;
        ylab = 'NL baseline param';
    case 'nl_slope'
        values = NL_params(:, 1)./(NL_params(:, 3)*sqrt(2/pi));
        ylims = [0 180];
        ytick = 0:50:150;
        ylab = 'NL slope';
    case 'area'
        values = rf_pixels*4.375^2; % in um^2
        ylims = [0 1.2e5];
        ytick = 0:6e4:1.2e5;
        ylab = 'RF area (\mum^2)';
    case 'diameter'
        values = avg_rad;
        ylims = [0 400];
        ytick = 0:200:400;
        ylab = 'RF diameter';
    case 'ellipse'
        values = elipse_ratio;
        ylims = [0 1];
        ytick = 0:0.5:1;
        ylab = 'RF ellipse ratio';
    case 'surround_center'
        values = surround_center;
        ylims = [0 0.05];
        ylab = 'Surround center ratio';
    case 'tftime2peak'
        values = -TF_time2peak;
        ylims = [0 120];
        ytick = 0:60:120;
        ylab = 'Time to peak(ms) abs';
    case 'tfwidth'
        values = TF_width;
        ylims = [30 70];
        ylab = 'Half width (ms)';
    case 'tfbiphasicpeaks'
        values = TF_biphasic_peaks;
        ylims = [0 0.9];
        ytick = 0:0.3:0.9;
        ylab = 'Biphasic index (peaks)';
    case 'tfbiphasicstregth'
        values = TF_biphasic_stregth;
        % ylims = [0 0.35];
        ylab = 'Biphasic index (strength)';

end
Davg = [];
Dsem = [];
clear xticlab
for i = 1:6
    Davg = [Davg mean(values(Ids{i}))];
    Dsem = [Dsem std(values(Ids{i}))/sqrt(num_n(i))];
    xticlab{i} = sprintf('%s (n=%d)', barlabels{i}, num_n(i));
end
box off
%%
p_ONn_t = ranksum(values(Ids{3}),values(Ids{4}))
p_ONFF_t = ranksum(values(Ids{5}),values(Ids{6}))
% [~, p_ONn_t] = ttest2(values(Ids{3}), values(Ids{4}))
% [~, p_ONFF_t] = ttest2(values(Ids{5}), values(Ids{6}))
%%
Colors = [0.3*ones(1, 3);
          0.6*ones(1, 3);
          [180 0 180]/255;
          [120  0 120]/255;
          [0  180 0]/255;
          [0 120 0]/255];
selection = 3:6;
figure; hold on
b = bar(selection, Davg(selection));
b.FaceColor = 'flat';
b.EdgeColor = 'w';
b.CData = Colors(selection, :);
errorbar(selection, Davg(selection), Dsem(selection), 'vertical', '|k', 'CapSize', 0');

% Overlay individual data points as dots
for i = 1:numel(selection)
    idx = selection(i);
    vals = values(Ids{idx});
    x_jitter = (rand(size(vals))-0.5)*0.15; % small horizontal jitter
    scatter(idx + x_jitter, vals, 40, 0.3*ones(1, 3), 'filled');
end

xticks(selection)
xticklabels(xticlab(selection))
clear y_ticklabels
if exist('ylims', 'var')
    if ~isempty(ylims)
        ylim(ylims);
        ylabel(ylab);
        yticks(ytick);
        for i = 1:length(ytick)
            y_ticklabels{i} = num2str(ytick(i));
        end
        yticklabels(y_ticklabels);
    end
end
xlim([2.5 6.5]);
%%
keyboard
%%
save_file_name = fullfile(save_folder, sprintf('BarPlot_RF%s_%s', eval_target, process_version(1:end-4)));
print(gcf, save_file_name, '-depsc', '-painters'); % EPS format
print(gcf, save_file_name, '-dpng', '-r300'); % PNG, 600 dpi

%%
Colors = parula(num_date);
eval_target = 'radius';
num_date = length(date_value);
date_id_array = nan(num_set, num_date);
for i = 1:num_date
    date_id_array(:, i) = date_ids == i;
end
num_n = sum(date_id_array, 1);
switch lower(eval_target)
    case 'radius'
        values = avg_rad;
        ylims = [50 120];
        ylab = 'RF radius';
    case 'ellipse'
        values = elipse_ratio;
        ylims = [0.6 1];
        ylab = 'RF ellipse ratio';
    case 'surround_center'
        values = surround_center;
        ylims = [0 0.05];
        ylab = 'Surround center ratio';
end
Davg = [];
Dsem = [];
clear xticlab
for i = 1:num_date
    Davg = [Davg mean(values(date_id_array(:, i)==1))];
    Dsem = [Dsem std(values(date_id_array(:, i)==1))/sqrt(num_n(i))];
    xticlab{i} = sprintf('%s (n=%d)', num2str(date_value(i)), num_n(i));
end
box off
figure; hold on
b = bar(1:num_date, Davg);
b.FaceColor = 'flat';
b.EdgeColor = 'w';
b.CData = Colors(1:num_date, :);
errorbar(1:num_date, Davg, Dsem, 'vertical', '|k', 'CapSize', 0');
xticks(1:num_date)
xticklabels(xticlab)
ylim(ylims);
ylabel(ylab);

