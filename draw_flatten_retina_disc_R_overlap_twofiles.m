clear; clc;
% close all;

% ================= USER SETTINGS =================
% Two Retistruct-style .mat files to overlay
% (must contain at least: x, y, phi0; and for rim_type='actual': phi, lambda, rim_idx)
dataFileA = 'C:\Users\jhsiang\Documents\R-retistruct\SMI32\data.mat';
dataFileB = 'C:\Users\jhsiang\Documents\R-retistruct\Mouse1R\data.mat';
figure_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';


% Colors for dataset A and B (RGB in [0..1])
colorA = [1  0 1];
colorB = [0  1 0];

% Rotations (degrees) used in the SECOND ROW only
rotation_ang_A = 0;
rotation_ang_B = 10;

% Toggle what to display: 'both' | 'only_A' | 'only_B'
showMode = 'both';

% Projection mode
method   = 'equidistant'; % 'equidistant' | 'equalarea' | 'conformal'
rim_type = 'actual';      % 'ideal' | 'actual'

% Density overlay appearance
gridRes = 220;
totalPrctile = 99;      % robust scale percentile for total density (A+B)
logKTotal = 0;       % log compression on total density (0 disables)
gammaTotal = 1.0;      % gamma on total brightness/value (<1 boosts background)
mixPower = 1;         % >1 tightens "equal mix" band around p=0.5
whiteStrength = 0;    % 0..1: how strongly overlap becomes white
whitePower = 1.7;       % >1 requires higher total density to whiten

% KDE bandwidth control
% Smaller -> more local (sharper) density; larger -> smoother density.
% Typical range: 0.3 to 1.5
kdeBandwidthScale = 0.7;

% Angle/scale overlays (coordinates are in degrees after rho_to_degrees)
showAngleScaleBar = true;
scaleBarDeg = 20;           % length of bar in degrees
showAngleGrid = true;      % rings/spokes on the disc
gridRingStepDeg = 30;       % concentric rings every N degrees
gridSpokeStepDeg = 180;      % spokes every N degrees
excludeSpokeAnglesDeg = [0 180];  % remove diagonal spokes (set [] to keep all)

% Slice settings (spherical space)
% Slice goes through the center along lambda = sliceAnglesDeg(1) and the opposite direction sliceAnglesDeg(2).
sliceAnglesDeg = [30 210];
nSliceSamples = 401;

% Slice stats: width at fraction of peak (e.g., 0.85 = drop to 85% of peak)
bandwidth_thr = 0.75;

% Unit conversion: mouse retina scale (visual degrees to retinal mm)
% If 1 visual degree corresponds to 32.5 um on the retina:
% KDE here outputs #/deg^2 (because xy are in degrees), so convert via:
%   (#/mm^2) = (#/deg^2) / (mmPerDeg^2)
micronsPerDeg = 32.5;
mmPerDeg = micronsPerDeg / 1000;
densityDeg2_to_mm2 = 1 / (mmPerDeg^2);

% ================= LOAD BOTH FILES =================
[A, nameA] = load_retistruct_file(dataFileA, rim_type);
[B, nameB] = load_retistruct_file(dataFileB, rim_type);

% Use dataset A as the common scaling reference (same circle)
phi0_common = A.phi0;

if isfield(B, 'phi0') && abs(B.phi0 - A.phi0) > 1e-9
    fprintf('Warning: phi0 differs between files (A=%.6g, B=%.6g). Using A as common scale.\n', A.phi0, B.phi0);
end

[projectionMode, isAreaPreserving, projectionLabel] = normalize_projection(method);

% Rim radius / extent based on common scale
rimRadius = phi_to_rho(phi0_common, phi0_common, isAreaPreserving);
plotExtent = 1.05 * rimRadius;
th  = linspace(0, 2*pi, 400);
idealCircle = rimRadius * [cos(th); sin(th)];

% ================= PROJECT (ROW 1: NO ROTATION) =================
[xyA_0, rimA_0] = project_one_dataset(A, phi0_common, projectionMode, isAreaPreserving, 0);
[xyB_0, rimB_0] = project_one_dataset(B, phi0_common, projectionMode, isAreaPreserving, 0);

% ================= PROJECT (ROW 2: WITH ROTATIONS) =================
[xyA_r, rimA_r] = project_one_dataset(A, phi0_common, projectionMode, isAreaPreserving, rotation_ang_A);
[xyB_r, rimB_r] = project_one_dataset(B, phi0_common, projectionMode, isAreaPreserving, rotation_ang_B);

% ================= KDE BANDWIDTH SUMMARY (CONSOLE) =================
% Coordinates are in visual degrees after rho_to_degrees().
% This prints the *base* bandwidth (per-axis) used in the pilot KDE and
% the distribution of *adaptive* per-point bandwidths.
print_kde_bandwidth_summary(xyA_0, kdeBandwidthScale, micronsPerDeg, sprintf('Dataset A (%s)', nameA));
print_kde_bandwidth_summary(xyB_0, kdeBandwidthScale, micronsPerDeg, sprintf('Dataset B (%s)', nameB));

% ================= FIGURE =================
figure('Color','w', 'Position', [100 100 1900 1000]);

% -------- Row 1 / Col 1: dots (overlaid) --------
subplot(2,3,1); hold on; axis equal;
set(gca, 'Visible', 'off');
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1);
plot(rimA_0(:,1), rimA_0(:,2), 'k-', 'LineWidth', 1.25);
if ~strcmpi(showMode, 'only_B')
    plot(xyA_0(:,1), xyA_0(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorA, 'MarkerEdgeColor', colorA);
end
if ~strcmpi(showMode, 'only_A')
    plot(xyB_0(:,1), xyB_0(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorB, 'MarkerEdgeColor', colorB);
end
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
add_angle_overlays(gca, rimRadius, plotExtent, showAngleScaleBar, scaleBarDeg, showAngleGrid, gridRingStepDeg, gridSpokeStepDeg, excludeSpokeAnglesDeg);
title(sprintf('Original overlay (%s)\nA=%s, B=%s', projectionLabel, nameA, nameB), 'Interpreter','none');

% -------- Row 1 / Col 2: density (overlaid RGB) --------
subplot(2,3,2); hold on; axis equal;
[gridX, gridY, insideMask] = make_disc_grid(plotExtent, rimRadius, gridRes);

dA0 = density_on_grid(xyA_0, gridX, gridY, insideMask, densityDeg2_to_mm2, kdeBandwidthScale);
dB0 = density_on_grid(xyB_0, gridX, gridY, insideMask, densityDeg2_to_mm2, kdeBandwidthScale);

switch lower(showMode)
    case 'only_a'
        dB0 = zeros(size(dB0));
    case 'only_b'
        dA0 = zeros(size(dA0));
    case 'both'
        % no-op
    otherwise
        error('Unknown showMode "%s". Use both | only_A | only_B.', showMode);
end

rgb0 = density_bivariate_to_rgb(dA0, dB0, insideMask, colorA, colorB, totalPrctile, logKTotal, gammaTotal, mixPower, whiteStrength, whitePower);
imagesc(gridX(1,:), gridY(:,1), rgb0);
set(gca, 'YDir', 'normal');
set(gca, 'Color', [1 1 1]);
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1.25);
plot(rimA_0(:,1), rimA_0(:,2), 'k-', 'LineWidth', 1.25);
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
axis off;
add_angle_overlays(gca, rimRadius, plotExtent, showAngleScaleBar, scaleBarDeg, showAngleGrid, gridRingStepDeg, gridSpokeStepDeg, excludeSpokeAnglesDeg);
title('Original overlay - density (bivariate: dominance hue, total brightness)');

% -------- Row 1 / Col 3: slice density profiles (A vs B) --------
subplot(2,3,3); hold on;
[sliceXY_0, sliceS_0] = make_spherical_diameter_slice(phi0_common, sliceAnglesDeg, projectionMode, isAreaPreserving, 0, nSliceSamples);
dSliceA0 = adaptiveKDE(xyA_0, sliceXY_0, kdeBandwidthScale);
dSliceB0 = adaptiveKDE(xyB_0, sliceXY_0, kdeBandwidthScale);

% adaptiveKDE returns a *probability density* (integrates to 1).
% Convert to *count density* by multiplying by number of cells, then convert units.
nCellsA0 = size(xyA_0, 1);
nCellsB0 = size(xyB_0, 1);
dSliceA0 = dSliceA0 * nCellsA0 * densityDeg2_to_mm2;
dSliceB0 = dSliceB0 * nCellsB0 * densityDeg2_to_mm2;

switch lower(showMode)
    case 'only_a'
        dSliceB0(:) = 0;
    case 'only_b'
        dSliceA0(:) = 0;
    case 'both'
        % no-op
end

plot(sliceS_0, dSliceA0, '-', 'Color', colorA, 'LineWidth', 2);
plot(sliceS_0, dSliceB0, '-', 'Color', colorB, 'LineWidth', 2);
grid off;
xlabel('Degrees from center (\x00B0)');
ylabel('KDE density (#/mm^2)');
maxAbs0 = max(abs(sliceS_0));
xlim([-maxAbs0, maxAbs0]);
ylim([0 60]);
yticks(0:30:60);
yticklabels({'0','30','60'});
legend({nameA, nameB}, 'Interpreter', 'none', 'Location', 'best');
title(sprintf('%d\x00B0 \x2194 %d\x00B0 slice density (original)', sliceAnglesDeg(1), sliceAnglesDeg(2)), 'Interpreter','none');

% -------- Row 2 / Col 1: dots (overlaid, rotated independently) --------
subplot(2,3,4); hold on; axis equal;
set(gca, 'Visible', 'off');
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1);
plot(rimA_r(:,1), rimA_r(:,2), 'k-', 'LineWidth', 1.25);
if ~strcmpi(showMode, 'only_B')
    plot(xyA_r(:,1), xyA_r(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorA, 'MarkerEdgeColor', colorA);
end
if ~strcmpi(showMode, 'only_A')
    plot(xyB_r(:,1), xyB_r(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorB, 'MarkerEdgeColor', colorB);
end
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
add_angle_overlays(gca, rimRadius, plotExtent, showAngleScaleBar, scaleBarDeg, showAngleGrid, gridRingStepDeg, gridSpokeStepDeg, excludeSpokeAnglesDeg);
title(sprintf('Rotated overlay (%s)\nA=%.1f\x00B0, B=%.1f\x00B0', projectionLabel, rotation_ang_A, rotation_ang_B), 'Interpreter','none');

% -------- Row 2 / Col 2: density (overlaid RGB, rotated independently) --------
subplot(2,3,5); hold on; axis equal;

dAr = density_on_grid(xyA_r, gridX, gridY, insideMask, densityDeg2_to_mm2, kdeBandwidthScale);
dBr = density_on_grid(xyB_r, gridX, gridY, insideMask, densityDeg2_to_mm2, kdeBandwidthScale);

switch lower(showMode)
    case 'only_a'
        dBr = zeros(size(dBr));
    case 'only_b'
        dAr = zeros(size(dAr));
    case 'both'
        % no-op
    otherwise
        error('Unknown showMode "%s". Use both | only_A | only_B.', showMode);
end

rgbr = density_bivariate_to_rgb(dAr, dBr, insideMask, colorA, colorB, totalPrctile, logKTotal, gammaTotal, mixPower, whiteStrength, whitePower);
imagesc(gridX(1,:), gridY(:,1), rgbr);
set(gca, 'YDir', 'normal');
set(gca, 'Color', [1 1 1]);
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1.25);
plot(rimA_r(:,1), rimA_r(:,2), 'k-', 'LineWidth', 1.25);
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
axis off;
add_angle_overlays(gca, rimRadius, plotExtent, showAngleScaleBar, scaleBarDeg, showAngleGrid, gridRingStepDeg, gridSpokeStepDeg, excludeSpokeAnglesDeg);
title('Rotated overlay - density (bivariate: dominance hue, total brightness)');

% -------- Row 2 / Col 3: slice density profiles (A vs B, rotated) --------
subplot(2,3,6); hold on;
[sliceXY_r, sliceS_r] = make_spherical_diameter_slice(phi0_common, sliceAnglesDeg, projectionMode, isAreaPreserving, rotation_ang_A, nSliceSamples);

% Evaluate on the same spherical slice for the rotated view:
% A uses rotation_ang_A, B uses rotation_ang_B
[sliceXY_Br, ~] = make_spherical_diameter_slice(phi0_common, sliceAnglesDeg, projectionMode, isAreaPreserving, rotation_ang_B, nSliceSamples);
dSliceAr = adaptiveKDE(xyA_r, sliceXY_r, kdeBandwidthScale);
dSliceBr = adaptiveKDE(xyB_r, sliceXY_Br, kdeBandwidthScale);

% Convert probability density -> count density, then deg^2 -> mm^2
nCellsAr = size(xyA_r, 1);
nCellsBr = size(xyB_r, 1);
dSliceAr = dSliceAr * nCellsAr * densityDeg2_to_mm2;
dSliceBr = dSliceBr * nCellsBr * densityDeg2_to_mm2;

switch lower(showMode)
    case 'only_a'
        dSliceBr(:) = 0;
    case 'only_b'
        dSliceAr(:) = 0;
    case 'both'
        % no-op
end

plot(sliceS_r, dSliceAr, '-', 'Color', colorA, 'LineWidth', 2);
plot(sliceS_r, dSliceBr, '-', 'Color', colorB, 'LineWidth', 2);
grid off;
xlabel('Degrees from center (\x00B0)');
ylabel('KDE density (#/mm^2)');
maxAbsr = max(abs(sliceS_r));
xlim([-maxAbsr, maxAbsr]);
ylim([0 60]);
yticks(0:30:60);
yticklabels({'0','30','60'});
legend({nameA, nameB}, 'Interpreter', 'none', 'Location', 'best');
title(sprintf('%d\x00B0 \x2194 %d\x00B0 slice density (rotated)', sliceAnglesDeg(1), sliceAnglesDeg(2)), 'Interpreter','none');

% --- Console stats for bottom-right slice curves ---
print_slice_curve_stats(sliceS_r, dSliceAr, dSliceBr, bandwidth_thr, micronsPerDeg, nameA, nameB, 'Bottom-right (rotated)');

% ================= SAVE FIGURE =================
% Create output folder if it doesn't exist
if ~exist(figure_folder, 'dir')
    mkdir(figure_folder);
    fprintf('Created folder: %s\n', figure_folder);
end

% Generate filename from current parameters
timestamp = datetime('now', 'Format', 'yyyyMMdd_HHmmss');
filename_base = sprintf('overlap_%s_rot%.0f_%s_gamma%.2f_mix%.1f_white%.1f_%s', ...
    showMode, rotation_ang_B, method, gammaTotal, mixPower, whiteStrength, timestamp);

% Save as PNG (high resolution, 300 dpi)
png_file = fullfile(figure_folder, [filename_base '.png']);
print(gcf, png_file, '-dpng', '-r300');
fprintf('Saved PNG: %s\n', png_file);

% Save as EPS (vector, Illustrator-compatible)
% Use -painters renderer for vector graphics and component editability
eps_file = fullfile(figure_folder, [filename_base '.eps']);
print(gcf, eps_file, '-depsc', '-painters', '-r300');
fprintf('Saved EPS: %s\n', eps_file);

fprintf('\nFiles saved to: %s\n', figure_folder);


% ================= HELPERS =================
function [Sout, shortName] = load_retistruct_file(dataFile, rim_type)
    S = load(dataFile);

    Sout.phi_data = S.x(:);
    Sout.lambda_data = S.y(:);

    if isfield(S, 'phi0')
        Sout.phi0 = S.phi0;
    else
        error('File %s missing phi0', dataFile);
    end

    switch rim_type
        case 'actual'
            if ~isfield(S, 'phi') || ~isfield(S, 'lambda') || ~isfield(S, 'rim_idx')
                error('For rim_type=actual, file %s must contain phi, lambda, rim_idx', dataFile);
            end
            Sout.phi_rim = S.phi(S.rim_idx);
            Sout.lambda_rim = S.lambda(S.rim_idx);
        case 'ideal'
            phi0  = Sout.phi0;
            theta = linspace(0, 2*pi, 400);
            Sout.phi_rim = phi0 * ones(size(theta));
            Sout.lambda_rim = theta - pi;
        otherwise
            error('Unknown rim_type "%s". Use "ideal" or "actual".', rim_type);
    end

    [~, base, ext] = fileparts(dataFile);
    shortName = [base ext];
end

function [data_xy, rim_xy] = project_one_dataset(S, phi0_common, projectionMode, isAreaPreserving, rotation_deg)
    phi_data = S.phi_data(:);
    lambda_data = S.lambda_data(:);
    phi_rim = S.phi_rim(:);
    lambda_rim = S.lambda_rim(:);

    if rotation_deg ~= 0
        rot_rad = rotation_deg * pi / 180;
        lambda_data = rotate_spherical_coords(lambda_data, rot_rad);
        lambda_rim  = rotate_spherical_coords(lambda_rim,  rot_rad);
    end

    [xr_raw, yr_raw] = sphere_spherical_to_polar_cart(phi_rim,  lambda_rim,  projectionMode);
    [xd_raw, yd_raw] = sphere_spherical_to_polar_cart(phi_data, lambda_data, projectionMode);

    rim_xy  = rho_to_degrees([xr_raw, yr_raw], phi0_common, isAreaPreserving);
    data_xy = rho_to_degrees([xd_raw, yd_raw], phi0_common, isAreaPreserving);
end

function [gridX, gridY, insideMask] = make_disc_grid(plotExtent, rimRadius, gridRes)
    gridLin = linspace(-plotExtent, plotExtent, gridRes);
    [gridX, gridY] = meshgrid(gridLin, gridLin);
    insideMask = hypot(gridX, gridY) <= (rimRadius + 1e-6);
end

function densityGrid = density_on_grid(points_xy, gridX, gridY, insideMask, densityDeg2_to_mm2, kdeBandwidthScale)
    densityGrid = nan(size(gridX));
    if any(insideMask(:))
        evalPoints = [gridX(insideMask), gridY(insideMask)];
        densityVals = adaptiveKDE(points_xy, evalPoints, kdeBandwidthScale);

        % adaptiveKDE returns a probability density (integrates to 1 over deg^2).
        % Convert to count density (#/mm^2).
        if nargin >= 5 && ~isempty(densityDeg2_to_mm2)
            nCells = size(points_xy, 1);
            densityVals = densityVals * nCells * densityDeg2_to_mm2;
        end
        densityGrid(insideMask) = densityVals;
    end
end

function rgb = density_bivariate_to_rgb(dA, dB, insideMask, colorA, colorB, totalPrctile, logKTotal, gammaTotal, mixPower, whiteStrength, whitePower)
% density_bivariate_to_rgb (Option B)
%   t = dA + dB controls brightness/value
%   p = dB / (dA + dB + eps) controls hue between colorA and colorB
%   If desired, strong overlap (p~0.5 AND high t) is pushed toward white.

    t = dA + dB;
    p = dB ./ (t + eps);

    % --- brightness from total density (robust scale + optional log compression) ---
    tvals = t(insideMask);
    tvals = tvals(~isnan(tvals));
    if isempty(tvals)
        s = 1;
    else
        s = prctile(tvals, totalPrctile);
        if ~isfinite(s) || s <= 0
            s = max(tvals);
        end
        if ~isfinite(s) || s <= 0
            s = 1;
        end
    end

    v = zeros(size(t));
    tin = max(0, t);
    if logKTotal > 0
        denom = log1p(logKTotal * s);
        if ~isfinite(denom) || denom <= 0
            denom = 1;
        end
        v(insideMask) = log1p(logKTotal * tin(insideMask)) ./ denom;
    else
        v(insideMask) = tin(insideMask) ./ max(s, eps);
    end
    v = min(1, max(0, v));
    v = v .^ gammaTotal;

    % --- hue from dominance ---
    base = zeros([size(t,1), size(t,2), 3]);
    for c = 1:3
        base(:,:,c) = (1 - p) .* colorA(c) + p .* colorB(c);
    end

    % --- whiten strong overlap only when both are present and total is high ---
    % mix = 1 at p=0.5 (equal), 0 at extremes
    mix = 1 - 2 * abs(p - 0.5);
    mix = min(1, max(0, mix));
    mix = mix .^ mixPower;

    w = whiteStrength .* (mix .* (v .^ whitePower));
    w = min(1, max(0, w));

    col = zeros([size(t,1), size(t,2), 3]);
    for c = 1:3
        col(:,:,c) = (1 - w) .* base(:,:,c) + w; % blend toward white
    end

    % --- render on white background with brightness v ---
    rgb = ones([size(t,1), size(t,2), 3]);
    for c = 1:3
        rgb(:,:,c) = 1 - v .* (1 - col(:,:,c));
    end

    for c = 1:3
        tmp = rgb(:,:,c);
        tmp(~insideMask) = 1;
        rgb(:,:,c) = tmp;
    end
end

function [x, y, rho] = sphere_spherical_to_polar_cart(phi, lambda, preserve)
    if nargin < 3 || isempty(preserve)
        preserve = 'latitude';
    end

    phi    = phi(:);
    lambda = lambda(:);

    preserve = lower(preserve);
    switch preserve
        case {'latitude', 'equidistant', 'lat'}
            rho = pi/2 + phi;
        case {'area', 'equalarea', 'pa'}
            rho = sqrt(2 * (1 + sin(phi)));
        case {'angle', 'conformal', 'stereographic'}
            rho = sqrt(2 * (1 + sin(phi)) ./ max(1 - sin(phi), eps));
        otherwise
            error('Unknown projection preserve mode "%s".', preserve);
    end

    x = rho .* cos(lambda);
    y = rho .* sin(lambda);
end

function lambda_rotated = rotate_spherical_coords(lambda, rotation_rad)
    lambda_rotated = lambda + rotation_rad;
    lambda_rotated = mod(lambda_rotated + pi, 2*pi) - pi;
end

function densityVals = adaptiveKDE(points, evalPoints, bandwidthScale)
    if nargin < 3 || isempty(bandwidthScale)
        bandwidthScale = 1.0;
    end
    bandwidthScale = max(bandwidthScale, eps);

    n = size(points, 1);
    if n < 2
        densityVals = zeros(size(evalPoints, 1), 1);
        return;
    end

    sigmaBase = std(points, 0, 1);
    sigmaBase(isnan(sigmaBase)) = 0;
    rangeVals = max(points, [], 1) - min(points, [], 1);
    zeroIdx = sigmaBase <= 0;
    sigmaBase(zeroIdx) = rangeVals(zeroIdx) / max(sqrt(n), 1);
    sigmaBase(sigmaBase <= 0) = max(rangeVals) * 0.01 + eps;

    baseBandwidth = 1.06 .* sigmaBase .* n^(-1/6);
    baseBandwidth = baseBandwidth .* bandwidthScale;
    baseBandwidth(baseBandwidth <= 0) = eps;

    pilot = mvksdensity(points, points, ...
                        'Bandwidth', baseBandwidth, ...
                        'Kernel', 'normal');
    pilot = max(pilot, realmin);
    geomMeanPilot = exp(mean(log(pilot)));
    lambda = (pilot / geomMeanPilot) .^ (-0.5);

    densityVals = zeros(size(evalPoints, 1), 1);
    for i = 1:n
        hi = baseBandwidth .* lambda(i);
        hi(hi <= 0) = eps;
        dx = (evalPoints(:, 1) - points(i, 1)) ./ hi(1);
        dy = (evalPoints(:, 2) - points(i, 2)) ./ hi(2);
        kernelVals = exp(-0.5 .* (dx.^2 + dy.^2)) ./ (2 * pi * hi(1) * hi(2));
        densityVals = densityVals + kernelVals;
    end

    densityVals = densityVals / n;
end

function print_slice_curve_stats(x, yA, yB, thrFrac, micronsPerDeg, nameA, nameB, label)
    if nargin < 4 || isempty(thrFrac)
        thrFrac = 0.85;
    end
    thrFrac = max(min(thrFrac, 1), 0);
    if nargin < 5
        micronsPerDeg = NaN;
    end
    if nargin < 6 || isempty(nameA)
        nameA = 'A';
    end
    if nargin < 7 || isempty(nameB)
        nameB = 'B';
    end
    if nargin < 8 || isempty(label)
        label = 'Slice stats';
    end

    statsA = peak_and_width_at_fraction(x, yA, thrFrac);
    statsB = peak_and_width_at_fraction(x, yB, thrFrac);

    % Cross-peak ratios
    ratioAatB = NaN;
    ratioBatA = NaN;
    if isfinite(statsA.peakY) && statsA.peakY > 0 && isfinite(statsB.peakX)
        yA_at_peakB = interp1_safe(x, yA, statsB.peakX);
        ratioAatB = yA_at_peakB / statsA.peakY;
    end
    if isfinite(statsB.peakY) && statsB.peakY > 0 && isfinite(statsA.peakX)
        yB_at_peakA = interp1_safe(x, yB, statsA.peakX);
        ratioBatA = yB_at_peakA / statsB.peakY;
    end

    fprintf('\n[SliceStats] %s (thr=%.2f of peak)\n', label, thrFrac);
    fprintf('  %s peak: x=%.3f deg, y=%.3f (#/mm^2)\n', nameA, statsA.peakX, statsA.peakY);
    fprintf('  %s width@thr (full): %.3f deg (L=%.3f, R=%.3f); half=%.3f deg\n', nameA, statsA.widthDeg, statsA.leftX, statsA.rightX, statsA.widthDeg/2);
    fprintf('  %s peak: x=%.3f deg, y=%.3f (#/mm^2)\n', nameB, statsB.peakX, statsB.peakY);
    fprintf('  %s width@thr (full): %.3f deg (L=%.3f, R=%.3f); half=%.3f deg\n', nameB, statsB.widthDeg, statsB.leftX, statsB.rightX, statsB.widthDeg/2);
    fprintf('  Cross-peak ratios: %s(x_peak_%s)/%s_peak = %.3f,  %s(x_peak_%s)/%s_peak = %.3f\n', ...
        nameA, nameB, nameA, ratioAatB, nameB, nameA, nameB, ratioBatA);

    if isfinite(micronsPerDeg)
        fprintf('  Width@thr (um): %s=%.1f, %s=%.1f\n', ...
            nameA, statsA.widthDeg * micronsPerDeg, nameB, statsB.widthDeg * micronsPerDeg);
    end
end

function s = peak_and_width_at_fraction(x, y, thrFrac)
    s = struct('peakX', NaN, 'peakY', NaN, 'leftX', NaN, 'rightX', NaN, 'widthDeg', NaN);
    if isempty(x) || isempty(y)
        return;
    end
    x = x(:);
    y = y(:);
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);
    if numel(x) < 3
        return;
    end

    % Ensure x is monotonic for well-defined left/right crossings.
    [x, ord] = sort(x);
    y = y(ord);

    [peakY, idxPeak] = max(y);
    if ~isfinite(peakY) || peakY <= 0
        s.peakX = x(idxPeak);
        s.peakY = peakY;
        return;
    end
    s.peakX = x(idxPeak);
    s.peakY = peakY;

    yThr = thrFrac * peakY;

    % Left crossing: find last index < idxPeak where y drops below threshold
    idxL = find(y(1:idxPeak) < yThr, 1, 'last');
    if ~isempty(idxL) && idxL < idxPeak
        s.leftX = interp_crossing_x(x(idxL), y(idxL), x(idxL+1), y(idxL+1), yThr);
    end

    % Right crossing: find first index > idxPeak where y drops below threshold
    idxRrel = find(y(idxPeak:end) < yThr, 1, 'first');
    if ~isempty(idxRrel)
        idxR = idxPeak + idxRrel - 1;
        if idxR > idxPeak
            s.rightX = interp_crossing_x(x(idxR-1), y(idxR-1), x(idxR), y(idxR), yThr);
        end
    end

    if isfinite(s.leftX) && isfinite(s.rightX)
        xL = min(s.leftX, s.rightX);
        xR = max(s.leftX, s.rightX);
        s.leftX = xL;
        s.rightX = xR;
        s.widthDeg = xR - xL;
    end
end

function xCross = interp_crossing_x(x1, y1, x2, y2, yThr)
    % Linear interpolation for x where y hits yThr between (x1,y1) and (x2,y2)
    if ~isfinite(x1) || ~isfinite(x2) || ~isfinite(y1) || ~isfinite(y2)
        xCross = NaN;
        return;
    end
    if y2 == y1
        xCross = (x1 + x2) / 2;
        return;
    end
    t = (yThr - y1) / (y2 - y1);
    xCross = x1 + t * (x2 - x1);
end

function yq = interp1_safe(x, y, xq)
    yq = NaN;
    if isempty(x) || isempty(y) || ~isfinite(xq)
        return;
    end
    x = x(:);
    y = y(:);
    good = isfinite(x) & isfinite(y);
    x = x(good);
    y = y(good);
    if numel(x) < 2
        return;
    end
    if xq < min(x) || xq > max(x)
        return;
    end
    yq = interp1(x, y, xq, 'linear');
end

function print_kde_bandwidth_summary(points, bandwidthScale, micronsPerDeg, label)
    if nargin < 2 || isempty(bandwidthScale)
        bandwidthScale = 1.0;
    end
    if nargin < 3 || isempty(micronsPerDeg)
        micronsPerDeg = NaN;
    end
    if nargin < 4 || isempty(label)
        label = 'Dataset';
    end

    n = size(points, 1);
    if n < 2
        fprintf('[KDE] %s: n=%d (not enough points for bandwidth estimate)\n', label, n);
        return;
    end

    sigmaBase = std(points, 0, 1);
    sigmaBase(isnan(sigmaBase)) = 0;
    rangeVals = max(points, [], 1) - min(points, [], 1);
    zeroIdx = sigmaBase <= 0;
    sigmaBase(zeroIdx) = rangeVals(zeroIdx) / max(sqrt(n), 1);
    sigmaBase(sigmaBase <= 0) = max(rangeVals) * 0.01 + eps;

    baseBW = 1.06 .* sigmaBase .* n^(-1/6);
    baseBW = baseBW .* max(bandwidthScale, eps);
    baseBW(baseBW <= 0) = eps;

    % Adaptive factors lambda(i) from pilot density at each data point
    pilot = mvksdensity(points, points, 'Bandwidth', baseBW, 'Kernel', 'normal');
    pilot = max(pilot, realmin);
    geomMeanPilot = exp(mean(log(pilot)));
    lambda = (pilot / geomMeanPilot) .^ (-0.5);

    hiX = baseBW(1) .* lambda;
    hiY = baseBW(2) .* lambda;

    p = [10 50 90];
    bx = prctile(hiX, p);
    by = prctile(hiY, p);

    fprintf('[KDE] %s: n=%d\n', label, n);
    fprintf('  Base BW (deg):  [%.3f  %.3f]\n', baseBW(1), baseBW(2));
    fprintf('  Adaptive BW (deg) 10/50/90%%:  x=[%.3f %.3f %.3f], y=[%.3f %.3f %.3f]\n', bx(1), bx(2), bx(3), by(1), by(2), by(3));

    if isfinite(micronsPerDeg)
        fprintf('  Base BW (um):   [%.1f  %.1f]\n', baseBW(1) * micronsPerDeg, baseBW(2) * micronsPerDeg);
        fprintf('  Adaptive BW (um) 10/50/90%%:  x=[%.1f %.1f %.1f], y=[%.1f %.1f %.1f]\n', ...
            bx(1) * micronsPerDeg, bx(2) * micronsPerDeg, bx(3) * micronsPerDeg, ...
            by(1) * micronsPerDeg, by(2) * micronsPerDeg, by(3) * micronsPerDeg);
    end
end

function pos = rho_to_degrees(pos, phi0, areaPreserving)
    if nargin < 3
        areaPreserving = false;
    end

    if areaPreserving
        phi0d = phi0 * 180 / pi;
        rho0  = spherical_to_polar_area(phi0);
        scale = (phi0d + 90) / rho0;
    else
        scale = 180 / pi;
    end

    pos = pos * scale;
end

function rho = phi_to_rho(phi, phi0, areaPreserving)
    if nargin < 3
        areaPreserving = false;
    end

    if areaPreserving
        phi0d = phi0 * 180 / pi;
        rho0  = spherical_to_polar_area(phi0);
        rho   = (phi0d + 90) / rho0 * spherical_to_polar_area(phi);
    else
        rho = (phi + pi/2) * 180 / pi;
    end
end

function rho = spherical_to_polar_area(phi)
    rho = sqrt(2 * (1 + sin(phi)));
end

function [mode, isAreaPreserving, label] = normalize_projection(method)
    if nargin < 1 || isempty(method)
        method = 'equidistant';
    end

    method = lower(method);
    switch method
        case {'equidistant', 'latitude', 'lat'}
            mode = 'latitude';
            isAreaPreserving = false;
            label = 'Azimuthal equidistant';
        case {'equalarea', 'area', 'pa'}
            mode = 'area';
            isAreaPreserving = true;
            label = 'Lambert azimuthal equal-area';
        case {'conformal', 'angle', 'stereographic'}
            mode = 'angle';
            isAreaPreserving = false;
            label = 'Stereographic (angle-preserving)';
        otherwise
            error('Unknown projection method "%s". Use equidistant | equalarea | conformal.', method);
    end
end

function add_angle_overlays(ax, rimRadiusDeg, plotExtent, showScaleBar, scaleBarDeg, showGrid, ringStepDeg, spokeStepDeg, excludeSpokeAnglesDeg)
    if nargin < 1 || isempty(ax) || ~isgraphics(ax)
        ax = gca;
    end

    hold(ax, 'on');

    % Optional polar grid (all in degrees)
    if showGrid
        % Concentric rings
        if nargin < 7 || isempty(ringStepDeg) || ringStepDeg <= 0
            ringStepDeg = 30;
        end
        ringR = ringStepDeg:ringStepDeg:floor(rimRadiusDeg / ringStepDeg) * ringStepDeg;
        th = linspace(0, 2*pi, 400);
        for r = ringR
            plot(ax, r*cos(th), r*sin(th), '-', 'Color', [0 0 0], 'LineWidth', 0.5);
            text(ax, 0, -r, sprintf('%g\x00B0', r), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 8, 'Color', [0 0 0]);
        end

        % Spokes
        if nargin < 8 || isempty(spokeStepDeg) || spokeStepDeg <= 0
            spokeStepDeg = 45;
        end
        if nargin < 9 || isempty(excludeSpokeAnglesDeg)
            excludeSpokeAnglesDeg = [];
        end
        spokes = 0:spokeStepDeg:(360 - spokeStepDeg);
        for ang = spokes
            if ~isempty(excludeSpokeAnglesDeg)
                angN = mod(ang, 360);
                exN = mod(excludeSpokeAnglesDeg(:)', 360);
                if any(abs(angN - exN) < 1e-9)
                    continue;
                end
            end
            a = deg2rad(ang);
            x = [0, rimRadiusDeg * cos(a)];
            y = [0, rimRadiusDeg * sin(a)];
            plot(ax, x, y, '-', 'Color', [0 0 0], 'LineWidth', 0.5);
        end
    end

    % Optional scale bar (in degrees)
    if showScaleBar
        if nargin < 5 || isempty(scaleBarDeg) || scaleBarDeg <= 0
            scaleBarDeg = 20;
        end

        y = -0.85 * rimRadiusDeg;
        x0 = -0.85 * rimRadiusDeg;
        x1 = x0 + scaleBarDeg;
        maxX = 0.85 * rimRadiusDeg;
        if x1 > maxX
            x0 = maxX - scaleBarDeg;
            x1 = maxX;
        end

        % Keep within plot limits even if rimRadiusDeg is close to plotExtent
        x0 = max(-plotExtent, min(plotExtent, x0));
        x1 = max(-plotExtent, min(plotExtent, x1));
        y  = max(-plotExtent, min(plotExtent, y));

        plot(ax, [x0 x1], [y y], 'k-', 'LineWidth', 2);
        plot(ax, [x0 x0], [y-0.8 y+0.8], 'k-', 'LineWidth', 1);
        plot(ax, [x1 x1], [y-0.8 y+0.8], 'k-', 'LineWidth', 1);
        text(ax, (x0+x1)/2, y, sprintf('%g\x00B0', scaleBarDeg), 'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 9, 'Color', 'k');
    end
end

function [sliceXY_deg, sliceS_deg] = make_spherical_diameter_slice(phi0, sliceAnglesDeg, projectionMode, isAreaPreserving, rotation_deg, nSamples)
    if nargin < 6 || isempty(nSamples)
        nSamples = 401;
    end
    if nargin < 5 || isempty(rotation_deg)
        rotation_deg = 0;
    end

    lam1 = deg2rad(sliceAnglesDeg(1));
    lam2 = deg2rad(sliceAnglesDeg(2));
    if rotation_deg ~= 0
        rot_rad = rotation_deg * pi / 180;
        lam1 = rotate_spherical_coords(lam1, rot_rad);
        lam2 = rotate_spherical_coords(lam2, rot_rad);
    end

    n1 = ceil(nSamples/2);
    n2 = nSamples - n1 + 1;

    % On this spherical cap, the center maps to phi = -pi/2 (rho=0).
    phi1 = linspace(phi0, -pi/2, n1)';
    phi2 = linspace(-pi/2, phi0, n2)';
    phi2 = phi2(2:end); % avoid duplicating the center sample

    phiSlice = [phi1; phi2];
    lamSlice = [repmat(lam1, n1, 1); repmat(lam2, numel(phi2), 1)];

    [xRaw, yRaw] = sphere_spherical_to_polar_cart(phiSlice, lamSlice, projectionMode);
    sliceXY_deg = rho_to_degrees([xRaw, yRaw], phi0, isAreaPreserving);

    % Signed position along the spherical diameter, in "degree" radial units.
    % This is rotation-invariant and matches the disc's radial scale.
    rhoDeg = phi_to_rho(phiSlice, phi0, isAreaPreserving);
    sliceS_deg = rhoDeg;
    sliceS_deg((n1+1):end) = -sliceS_deg((n1+1):end);
end
