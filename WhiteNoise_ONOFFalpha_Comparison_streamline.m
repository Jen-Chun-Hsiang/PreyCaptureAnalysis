% WhiteNoise_ONOFFalpha_Comparison_streamline.m
%
% Streamlined version of WhiteNoise_ONOFFalpha_Comparison.m
% - Optionally loads an existing processed .mat via file dialog (uigetfile)
% - If no processed file is selected/found, runs TF fitting and saves processed output
% - Centralizes figure generation and statistical tests
%
% This script intentionally preserves the core computations (TF metrics and
% group definitions) used in WhiteNoise_ONOFFalpha_Comparison.m.

close all; clear; clc;

% Ensure deterministic plots/stat outputs across runs (e.g., jittered scatter points)
rng(0, 'twister');

%% ----------------------- User configuration -----------------------
cfg = struct();

% Data + metadata
cfg.data_sets = {'e100724', 'f100724', 'a101224', 'b101224',  'c101224',  'd101424', 'e101224',...
                 'b101424', 'c101424', 'd101424', 'e101424',  'a101624',  'b101624', 'd101624', 'e101624',...
                 'b101924', 'c101924', 'd101924', 'e101924',  'b103124',  'e103124', 'a110424', 'b110424',...
                 'c110424', 'd110424', 'e110424', 'f110424',  'g110424',  'a110924', 'b110924', 'c110924',...
                 'a111224'};

cfg.cell_type = {'OFF',      'OFF',    'OFF',      'ON',       'OFF',     'ON',      'OFF',...
                 'OFF',      'OFF',    'ON',       'ON',       'ON',      'ON',      'ON',       'ON',...
                 'ON',       'OFF',    'OFF',      'OFF',      'ON',      'OFF',     'ON',       'ON',...
                 'ON',       'ON',     'OFF',      'ON',       'OFF',     'ON',      'OFF',      'OFF',...
                 'ON'};

cfg.location =  {'Temporal', 'Temporal','Nasal',   'Nasal',    'Nasal',   'Nasal',   'Nasal',...
                 'Temporal', 'Temporal','Temporal','Temporal', 'Nasal',   'Nasal',   'Nasal',    'Nasal',...
                 'Nasal',    'Nasal',   'Nasal',   'Nasal',    'Temporal','Temporal','Temporal', 'Temporal',...
                 'Temporal', 'Temporal','Temporal','Temporal', 'Temporal','Temporal','Temporal', 'Temporal',...
                 'Temporal'};

% Where per-cell .mat files live (stdSTA, tRF)
cfg.folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';

% Processing / fit choices
cfg.is_normalized_tf = true;       % normalize csig by max(abs(csig)) before fitting/plotting
cfg.temporalFitFcn = @GaussianTemporalFilter2; % keep consistent with WhiteNoise_ONOFFalpha_Comparison.m

% Time base
cfg.Fz = 100;
cfg.WinT = [-0.5 0];

% I/O
cfg.save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';
if ~exist(cfg.save_folder, 'dir'), mkdir(cfg.save_folder); end
cfg.save_tf_folder = fullfile(cfg.save_folder, 'TF_Fits');
if ~exist(cfg.save_tf_folder, 'dir'), mkdir(cfg.save_tf_folder); end

% Processed file options
cfg.process_version_default = 'GaussianFitting_processed_082025_1.mat';
cfg.use_file_dialog_to_load_processed = true;   % if true -> uigetfile prompt
cfg.auto_save_processed_if_computed = true;
% If you LOAD an existing processed file, keep it read-only by default.
% (Set true only if you intentionally want this script to modify the selected .mat.)
cfg.allow_append_missing_fields_to_loaded_processed = false;

% Plot + stats targets
% Choose one: 'area' | 'diameter' | 'tftime2peak' | 'tfwidth' | 'tfbiphasicpeaks' | 'tfbiphasicstregth'
% plus RF geometry vars if present in processed file.
cfg.eval_target = 'area';

% Stats configuration
cfg.alpha = 0.05;

% RF geometry constants (match original script assumptions)
cfg.pixel_um = 4.375;  % microns per pixel
cfg.rf_threshold_std_plus = 2; % used only if rf_pixels must be computed

%% ----------------------- Select + load processed output early -----------------------
if cfg.use_file_dialog_to_load_processed
    [f, p] = uigetfile({'*.mat','MAT-files (*.mat)'}, 'Select processed GaussianFitting_processed_*.mat', cfg.folder_name);
    if isequal(f, 0)
        % Fall back to default
        processedFile = fullfile(cfg.folder_name, cfg.process_version_default);
    else
        processedFile = fullfile(p, f);
    end
else
    processedFile = fullfile(cfg.folder_name, cfg.process_version_default);
end

hasProcessed = exist(processedFile, 'file') == 2;
S = struct();
if hasProcessed
    disp(['Loading processed file: ' processedFile]);
    S = load(processedFile);
end

%% ----------------------- Load per-cell data -----------------------
Data = cell(numel(cfg.data_sets), 1);
for i = 1:numel(cfg.data_sets)
    file_name = sprintf('%s.mat', cfg.data_sets{i});
    Data{i} = load(fullfile(cfg.folder_name, file_name), 'stdSTA', 'tRF');
end

% Convert types/locations to numeric
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), cfg.location);
cell_type_numeric = cellfun(@(x) strcmp(x, 'ON'), cfg.cell_type);

% Common time axis for temporal filter trace
t = cfg.WinT(1):1/cfg.Fz:cfg.WinT(end);
ct = t(2:end);

%% ----------------------- Load or compute processed output -----------------------
if ~hasProcessed
    disp('No processed file found/selected. Computing temporal fits now...');
    S = struct();

    % Fit temporal filters
    Trace = nan(numel(cfg.data_sets), numel(ct));
    x = 1:numel(ct);

    for k = 1:numel(cfg.data_sets)
        fprintf('TF fitting %s... %d/%d\n', cfg.data_sets{k}, k, numel(cfg.data_sets));
        csig = Data{k}.tRF;
        if cfg.is_normalized_tf
            csig = csig ./ max(abs(csig));
        end

        OptW = cfg.temporalFitFcn(csig');

        tf = gaussmf(x, [OptW(1) OptW(3)])*OptW(5) ...
           - gaussmf(x, [OptW(2) OptW(4)])*(OptW(6)*OptW(5)) + OptW(7);

        if k == 1
            Gauss_TF_est = nan(numel(cfg.data_sets), numel(OptW));
        end
        Gauss_TF_est(k, :) = OptW;
        Trace(k, :) = csig;

        % Save a per-cell fit figure (optional, but keeps parity with original script)
        fig = figure('Visible','off'); hold on
        plot(ct, csig, 'k', 'LineWidth', 1.5);
        plot(ct, tf, 'r--', 'LineWidth', 1.5);
        legend('Normalized csig', 'Fitted TF');
        xlabel('Time (s)'); ylabel('Normalized Value');
        title(sprintf('TF Fit: %s (%s, %s)', cfg.data_sets{k}, cfg.cell_type{k}, cfg.location{k}));
        set(gca, 'Box', 'off');
        saveas(fig, fullfile(cfg.save_tf_folder, sprintf('TFfit_%s.png', cfg.data_sets{k})));
        close(fig);
    end

    S.ct = ct;
    S.Trace = Trace;
    S.Gauss_TF_est = Gauss_TF_est;
    S.location_type_numeric = location_type_numeric;
    S.cell_type_numeric = cell_type_numeric;
    S.cell_type = cfg.cell_type;
    S.location = cfg.location;
    S.data_sets = cfg.data_sets;

    if cfg.auto_save_processed_if_computed
        disp(['Saving processed file: ' processedFile]);
        save(processedFile, '-struct', 'S');
    end
end

% Normalize and validate needed fields
if ~isfield(S, 'ct'), S.ct = ct; end
if ~isfield(S, 'location_type_numeric'), S.location_type_numeric = location_type_numeric; end
if ~isfield(S, 'cell_type_numeric'), S.cell_type_numeric = cell_type_numeric; end
if ~isfield(S, 'data_sets'), S.data_sets = cfg.data_sets; end

%% ----------------------- Compute temporal metrics (matches original) -----------------------
% This uses Data{i}.tRF (raw) + S.Gauss_TF_est (fit params)
TF_time2peak = nan(numel(cfg.data_sets), 1);
TF_width = nan(numel(cfg.data_sets), 1);
TF_biphasic_peaks = nan(numel(cfg.data_sets), 1);
TF_biphasic_stregth = nan(numel(cfg.data_sets), 1);

hwith_thr = 0.5;
nt = size(S.Trace, 2);
isON = S.cell_type_numeric == 1;

for i = 1:numel(cfg.data_sets)
    csig = Data{i}.tRF(:)';

    % Interpolate for finer resolution (as in original)
    csig_interp = interp1(linspace(0, 1, nt), csig, linspace(0, 1, 1000), 'cubic');
    t_interp = linspace(cfg.WinT(1), cfg.WinT(end), 1000);

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
        TF_biphasic_stregth(i) = S.Gauss_TF_est(i, 6)./S.Gauss_TF_est(i, 5);
    else
        csig_thr = csig_interp < -hwith_thr;
        TF_biphasic_peaks(i) = 1-2*abs(minv/(maxv+minv)-0.5);
        TF_biphasic_stregth(i) = S.Gauss_TF_est(i, 5)./S.Gauss_TF_est(i, 6);
    end

    TF_width(i) = sum(csig_thr)*(1/cfg.Fz)*nt;
end

%% ----------------------- Plot: ON/OFF average temporal filters -----------------------
figure('Name','Temporal filters (group averages)'); hold on
TracePlot = S.Trace;
if cfg.is_normalized_tf
    yt = -1:0.5:1;
else
    yt = -100:50:100;
end

h1 = plot(S.ct, squeeze(mean(TracePlot(S.cell_type_numeric==1, :), 1)), 'Color', [245 182 66]/255, 'LineWidth', 2);
h2 = plot(S.ct, squeeze(mean(TracePlot(S.cell_type_numeric==0, :), 1)), 'Color', 0*ones(1, 3), 'LineWidth', 2);
xlabel('Time (s)'); ylabel('Average stimulus value');
xticks(-0.5:0.25:0);
xticklabels({'-0.5','-0.25','0'});
yticks(yt);
legend([h1,h2], 'ON', 'OFF');
set(gca,'Box','off');

%% ----------------------- Plot: Biphasy colored traces + location means -----------------------
Disp_Type = 'ON';
cell_type_id = strcmpi(Disp_Type, 'ON');
type_ids = find(S.cell_type_numeric == cell_type_id)';

biphase_for_color = S.Gauss_TF_est(type_ids, 5:6);
if cell_type_id == 1
    biphase_for_color = biphase_for_color(:, 2)./biphase_for_color(:, 1);
else
    biphase_for_color = biphase_for_color(:, 1)./biphase_for_color(:, 2);
end

colors = parula(256);
color_ids = round(((biphase_for_color-min(biphase_for_color))*255/range(biphase_for_color))+1);
color_ids(isnan(color_ids)) = 1;

fig = figure('Name', sprintf('TF traces colored by biphasy (%s)', Disp_Type));
subplot(1, 2, 1); hold on
for k = 1:numel(type_ids)
    csig = Data{type_ids(k)}.tRF;
    if cfg.is_normalized_tf
        csig = csig./max(abs(csig));
    end
    plot(S.ct, csig(:), 'Color', colors(color_ids(k), :));
end
xlabel('Time (s)'); ylabel('Normalized TF'); set(gca,'Box','off'); title('Traces (colored by biphasy)');

subplot(1, 2, 2); hold on
for k = 1:numel(type_ids)
    cid = type_ids(k);
    csig = Data{cid}.tRF;
    if cfg.is_normalized_tf
        csig = csig./max(abs(csig));
    end
    if strcmpi(S.location{cid}, 'temporal')
        plot(S.ct, csig(:), 'Color', [66 182 245]/255);
    else
        plot(S.ct, csig(:), 'Color', [247 153 205]/255);
    end
end
hT = plot(S.ct, squeeze(mean(TracePlot(S.cell_type_numeric==cell_type_id & S.location_type_numeric == 1, :), 1)),...
    'Color', [27 59 242]/255, 'LineWidth', 2);
hN = plot(S.ct, squeeze(mean(TracePlot(S.cell_type_numeric==cell_type_id & S.location_type_numeric == 0, :), 1)),...
    'Color', [242 27 145]/255, 'LineWidth', 2);
legend([hT,hN], 'Temporal', 'Nasal');
xlabel('Time (s)'); ylabel('Normalized TF'); set(gca,'Box','off'); title('By location');

%% ----------------------- Stats: location comparison within Disp_Type -----------------------
% NOTE: biphase_for_color is a *fit-based ratio* derived from S.Gauss_TF_est.
% To avoid confusion, only report this location comparison when the user-selected
% cfg.eval_target matches this quantity.
if strcmpi(cfg.eval_target, 'tfbiphasicstregth')
    locs = S.location_type_numeric(type_ids);
    x_nasal = biphase_for_color(locs==0);
    y_temporal = biphase_for_color(locs==1);

    p_biphasy_loc_ttest = NaN;
    if numel(x_nasal) >= 2 && numel(y_temporal) >= 2
        [~, p_biphasy_loc_ttest] = ttest2(x_nasal, y_temporal);
    end

    fprintf('\n[%s] %s Nasal vs Temporal: ttest2 p=%.4g\n', ...
        cfg.eval_target, Disp_Type, p_biphasy_loc_ttest);
end

%% ----------------------- Bar plot + stats for selected eval_target -----------------------
Ids = cell(6,1);
Ids{1} = S.cell_type_numeric == 1;
Ids{2} = S.cell_type_numeric == 0;
Ids{3} = S.cell_type_numeric == 1 & S.location_type_numeric == 1;
Ids{4} = S.cell_type_numeric == 1 & S.location_type_numeric == 0;
Ids{5} = S.cell_type_numeric == 0 & S.location_type_numeric == 1;
Ids{6} = S.cell_type_numeric == 0 & S.location_type_numeric == 0;
barlabels = {'ON', 'OFF', 'ON-temporal', 'ON-nasal', 'OFF-temporal', 'OFF-nasal'};
num_n = cellfun(@sum, Ids);

% ---- RF geometry metrics (loaded from processed file; computed only if needed) ----
gauss_est = [];
rf_pixels = [];
avg_rad = [];
elipse_ratio = [];
surround_center = [];

if isfield(S, 'gauss_est')
    gauss_est = S.gauss_est;
    if size(gauss_est, 2) >= 4
        elipse_ratio = min(gauss_est(:, 3:4), [], 2) ./ max(gauss_est(:, 3:4), [], 2);
        avg_rad = 2*sqrt(gauss_est(:, 3).*gauss_est(:, 4)) * 2 * cfg.pixel_um;
    end
    if size(gauss_est, 2) >= 10
        surround_center = abs(gauss_est(:, 10) ./ gauss_est(:, 7));
    end
end

if isfield(S, 'rf_pixels')
    rf_pixels = S.rf_pixels;
end

% Enforce required fields for RF-geometry targets (avoid silent all-NaN outputs)
if any(strcmpi(cfg.eval_target, {'diameter','ellipse','surround_center'})) && isempty(gauss_est)
    error(['eval_target="' cfg.eval_target '" requires gauss_est in the processed file. ' ...
           'Load a GaussianFitting_processed_*.mat generated by WhiteNoise_ONOFFalpha_Comparison.m (with spatial fit).']);
end

% If user requests RF area but rf_pixels isn't in the processed file, compute it quickly
% from gauss_est + stdSTA (matches FindThreshold_MeanMinusKStd.m behavior).
if strcmpi(cfg.eval_target, 'area') && isempty(rf_pixels)
    if isempty(gauss_est)
        error(['eval_target="area" requires rf_pixels or gauss_est in the processed file. ' ...
               'Load a GaussianFitting_processed_*.mat generated by WhiteNoise_ONOFFalpha_Comparison.m (with spatial fit), ' ...
               'or re-run that script once to generate gauss_est/rf_pixels.']);
    end
    rf_pixels = compute_rf_pixels_from_gauss(Data, gauss_est, cfg.rf_threshold_std_plus);

    % Only append to disk if explicitly allowed.
    if (~hasProcessed) || cfg.allow_append_missing_fields_to_loaded_processed
        try
            save(processedFile, 'rf_pixels', '-append');
            disp('Appended rf_pixels to processed file for future runs.');
        catch me
            warning('WhiteNoise:AppendProcessedFailed', '%s', sprintf('Could not append rf_pixels to processed file: %s', me.message));
        end
    else
        disp('Computed rf_pixels in-memory (not appended to loaded processed file).');
    end
end

[values, ylab, ylims, ytick] = resolve_eval_target(cfg.eval_target, ...
    TF_time2peak, TF_width, TF_biphasic_peaks, TF_biphasic_stregth, ...
    rf_pixels, avg_rad, elipse_ratio, surround_center, cfg.pixel_um);

%% ----------------------- Console summary stats -----------------------
fprintf('\n[%s] Summary stats by group (mean, std)\n', cfg.eval_target);
for gi = 1:6
    v = values(Ids{gi});
    v = v(~isnan(v));
    n = numel(v);
    if n == 0
        fprintf('  %s: n=0 (all NaN or empty)\n', barlabels{gi});
    else
        fprintf('  %s: n=%d, mean=%.4g, std=%.4g\n', barlabels{gi}, n, mean(v), std(v));
    end
end

Davg = nan(1,6);
Dsem = nan(1,6);
xticlab = cell(1,6);
for i = 1:6
    % Match original script aggregation exactly
    v = values(Ids{i});
    Davg(i) = mean(v);
    Dsem(i) = std(v) / sqrt(max(1, num_n(i)));
    xticlab{i} = sprintf('%s (n=%d)', barlabels{i}, num_n(i));
end

figure('Name', ['Bar plot - ' cfg.eval_target]); hold on
selection = 3:6;
Colors = [0.3*ones(1, 3);
          0.6*ones(1, 3);
          [180 0 180]/255;
          [120  0 120]/255;
          [0  180 0]/255;
          [0 120 0]/255];

b = bar(selection, Davg(selection));
b.FaceColor = 'flat';
b.EdgeColor = 'w';
b.CData = Colors(selection, :);
errorbar(selection, Davg(selection), Dsem(selection), 'vertical', '|k', 'CapSize', 0');

for i = 1:numel(selection)
    idx = selection(i);
    vals = values(Ids{idx});
    x_jitter = (rand(size(vals))-0.5)*0.15;
    scatter(idx + x_jitter, vals, 40, 0.3*ones(1, 3), 'filled');
end

xticks(selection);
xticklabels(xticlab(selection));
xlim([2.5 6.5]);
if ~isempty(ylims), ylim(ylims); end
if ~isempty(ylab), ylabel(ylab); end
if ~isempty(ytick), yticks(ytick); end
set(gca,'Box','off');

% Stats: ON vs OFF (ignore location)
v_on = values(Ids{1}); v_on = v_on(~isnan(v_on));
v_off = values(Ids{2}); v_off = v_off(~isnan(v_off));
p_onoff_ttest = NaN;
if numel(v_on) >= 2 && numel(v_off) >= 2
    [~, p_onoff_ttest] = ttest2(v_on, v_off);
end
fprintf('\n[%s] ON vs OFF: ttest2 p=%.4g\n', cfg.eval_target, p_onoff_ttest);

% Stats: Nasal vs Temporal (ignore ON/OFF)
v_nasal_all = values(S.location_type_numeric == 0); v_nasal_all = v_nasal_all(~isnan(v_nasal_all));
v_temporal_all = values(S.location_type_numeric == 1); v_temporal_all = v_temporal_all(~isnan(v_temporal_all));
p_nt_ttest = NaN;
if numel(v_nasal_all) >= 2 && numel(v_temporal_all) >= 2
    [~, p_nt_ttest] = ttest2(v_nasal_all, v_temporal_all);
end
fprintf('[%s] Nasal vs Temporal (all cells): ttest2 p=%.4g\n', cfg.eval_target, p_nt_ttest);

% Stats between temporal vs nasal within ON and within OFF
fprintf('[%s] Stats (ttest2): ON temporal vs ON nasal\n', cfg.eval_target);
[~, p_ONn_t_ttest] = ttest2(values(Ids{3}), values(Ids{4}));
disp(p_ONn_t_ttest);

fprintf('[%s] Stats (ttest2): OFF temporal vs OFF nasal\n', cfg.eval_target);
[~, p_OFFn_t_ttest] = ttest2(values(Ids{5}), values(Ids{6}));
disp(p_OFFn_t_ttest);

%% ----------------------- Local helper functions -----------------------
function [values, ylab, ylims, ytick] = resolve_eval_target(eval_target, TF_time2peak, TF_width, TF_biphasic_peaks, TF_biphasic_stregth, rf_pixels, avg_rad, elipse_ratio, surround_center, pixel_um)
    ylims = []; ytick = []; ylab = ''; %#ok<NASGU>
    switch lower(eval_target)
        case 'area'
            values = rf_pixels * pixel_um^2; % um^2
            ylims = [0 1.5e5];
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
            values = -TF_time2peak; % keep sign convention from original script
            ylims = [0 120];
            ytick = 0:60:120;
            ylab = 'Time to peak (ms) abs';
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
            ylab = 'Biphasic index (strength)';
        otherwise
            error('Unknown eval_target: %s', eval_target);
    end
end

function rf_pixels = compute_rf_pixels_from_gauss(Data, gauss_est, std_plus)
    % Replicates FindThreshold_MeanMinusKStd.m (but as a function; no workspace dependence).
    num_set = size(gauss_est, 1);
    rf_pixels = nan(num_set, 1);
    for k = 1:num_set
        image = Data{k}.stdSTA';
        params = gauss_est(k, :);
        cx = params(1); cy = params(2); sx = params(3); sy = params(4); theta = params(5);

        [X, Y] = meshgrid(1:size(image, 2), 1:size(image, 1));
        Xc = X - cx; Yc = Y - cy;
        Xr = Xc * cos(theta) + Yc * sin(theta);
        Yr = -Xc * sin(theta) + Yc * cos(theta);
        mask = (Xr.^2/(2*sx^2) + Yr.^2/(2*sy^2)) <= 2^2; % 2*sigma ellipse

        mask_pixels = image(~mask);
        mu = mean(mask_pixels);
        sigma = std(mask_pixels);
        threshold = mu + std_plus*sigma;

        bw = image > threshold;
        CC = bwconncomp(bw);
        stats = regionprops(CC, 'Area');
        if isempty(stats)
            largest_area = 0;
        else
            areas = [stats.Area];
            largest_area = max(areas);
        end
        rf_pixels(k) = largest_area;
    end
end
