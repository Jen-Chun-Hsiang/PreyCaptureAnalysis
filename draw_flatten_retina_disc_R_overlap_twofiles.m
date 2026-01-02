clear; clc;
% close all;

% ================= USER SETTINGS =================
% Two Retistruct-style .mat files to overlay
% (must contain at least: x, y, phi0; and for rim_type='actual': phi, lambda, rim_idx)
dataFileA = 'C:\Users\jhsiang\Documents\R-retistruct\SMI32\data.mat';
dataFileB = 'C:\Users\jhsiang\Documents\R-retistruct\Mouse1R\data.mat';

% Colors for dataset A and B (RGB in [0..1])
colorA = [1 1 0];
colorB = [0 0 1];

% Rotations (degrees) used in the SECOND ROW only
rotation_ang_A = 0;
rotation_ang_B = 10;

% Projection mode
method   = 'equidistant'; % 'equidistant' | 'equalarea' | 'conformal'
rim_type = 'actual';      % 'ideal' | 'actual'

% Density overlay appearance
gridRes = 220;
gammaA  = 0.65;
gammaB  = 0.65;
alphaMax = 0.95; % overall alpha cap for density overlay

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

% ================= FIGURE =================
figure('Color','w', 'Position', [100 100 1400 1000]);

% -------- Row 1 / Col 1: dots (overlaid) --------
subplot(2,2,1); hold on; axis equal;
set(gca, 'Visible', 'off');
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1);
plot(rimA_0(:,1), rimA_0(:,2), 'k-', 'LineWidth', 1.25);
plot(xyA_0(:,1), xyA_0(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorA, 'MarkerEdgeColor', colorA);
plot(xyB_0(:,1), xyB_0(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorB, 'MarkerEdgeColor', colorB);
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
title(sprintf('Original overlay (%s)\nA=%s, B=%s', projectionLabel, nameA, nameB), 'Interpreter','none');

% -------- Row 1 / Col 2: density (overlaid RGB) --------
subplot(2,2,2); hold on; axis equal;
[gridX, gridY, insideMask] = make_disc_grid(plotExtent, rimRadius, gridRes);

dA0 = density_on_grid(xyA_0, gridX, gridY, insideMask);
dB0 = density_on_grid(xyB_0, gridX, gridY, insideMask);

rgb0 = density_pair_to_rgb(dA0, dB0, insideMask, colorA, colorB, gammaA, gammaB);
alpha0 = min(alphaMax, 0.15 + 0.85 * normalize01(nansum(cat(3,dA0,dB0),3), insideMask));

imagesc(gridX(1,:), gridY(:,1), rgb0, 'AlphaData', alpha0);
set(gca, 'YDir', 'normal');
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1.25);
plot(rimA_0(:,1), rimA_0(:,2), 'k-', 'LineWidth', 1.25);
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
axis off;
title('Original overlay - density (RGB blend)');

% -------- Row 2 / Col 1: dots (overlaid, rotated independently) --------
subplot(2,2,3); hold on; axis equal;
set(gca, 'Visible', 'off');
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1);
plot(rimA_r(:,1), rimA_r(:,2), 'k-', 'LineWidth', 1.25);
plot(xyA_r(:,1), xyA_r(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorA, 'MarkerEdgeColor', colorA);
plot(xyB_r(:,1), xyB_r(:,2), 'o', 'MarkerSize', 3, 'MarkerFaceColor', colorB, 'MarkerEdgeColor', colorB);
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
title(sprintf('Rotated overlay (%s)\nA=%.1f\x00B0, B=%.1f\x00B0', projectionLabel, rotation_ang_A, rotation_ang_B), 'Interpreter','none');

% -------- Row 2 / Col 2: density (overlaid RGB, rotated independently) --------
subplot(2,2,4); hold on; axis equal;

dAr = density_on_grid(xyA_r, gridX, gridY, insideMask);
dBr = density_on_grid(xyB_r, gridX, gridY, insideMask);

rgbr = density_pair_to_rgb(dAr, dBr, insideMask, colorA, colorB, gammaA, gammaB);
alphar = min(alphaMax, 0.15 + 0.85 * normalize01(nansum(cat(3,dAr,dBr),3), insideMask));

imagesc(gridX(1,:), gridY(:,1), rgbr, 'AlphaData', alphar);
set(gca, 'YDir', 'normal');
plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1.25);
plot(rimA_r(:,1), rimA_r(:,2), 'k-', 'LineWidth', 1.25);
xlim([-plotExtent, plotExtent]);
ylim([-plotExtent, plotExtent]);
axis off;
title('Rotated overlay - density (RGB blend)');


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

function densityGrid = density_on_grid(points_xy, gridX, gridY, insideMask)
    densityGrid = nan(size(gridX));
    if any(insideMask(:))
        evalPoints = [gridX(insideMask), gridY(insideMask)];
        densityVals = adaptiveKDE(points_xy, evalPoints);
        densityGrid(insideMask) = densityVals;
    end
end

function rgb = density_pair_to_rgb(dA, dB, insideMask, colorA, colorB, gammaA, gammaB)
    a = normalize01(dA, insideMask) .^ gammaA;
    b = normalize01(dB, insideMask) .^ gammaB;

    % Multiply-blend on white background (keeps both colors visible)
    % rgb = 1 - (1 - colorA*a) * (1 - colorB*b)
    rgb = ones([size(dA,1), size(dA,2), 3]);
    for c = 1:3
        layerA = 1 - colorA(c) .* a;
        layerB = 1 - colorB(c) .* b;
        rgb(:,:,c) = 1 - (layerA .* layerB);
    end

    % Outside disc: pure white
    for c = 1:3
        tmp = rgb(:,:,c);
        tmp(~insideMask) = 1;
        rgb(:,:,c) = tmp;
    end
end

function out = normalize01(v, insideMask)
    out = zeros(size(v));
    mx = max(v(insideMask), [], 'omitnan');
    if isempty(mx) || ~isfinite(mx) || mx <= 0
        mx = 1;
    end
    out(insideMask) = v(insideMask) ./ mx;
    out(~insideMask) = 0;
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

function densityVals = adaptiveKDE(points, evalPoints)
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
