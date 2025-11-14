clear; close all; clc;
dataFile = 'C:\Users\jhsiang\Documents\R-retistruct\Test07\data.mat';



% DIAGNOSTIC: Inspect data ranges and distribution
S = load(dataFile);
phi_data    = S.x;
lambda_data = S.y;
rim_type = 'actual'; % 'ideal' | 'actual'

switch rim_type
    case 'actual'
        phi0 = S.phi0;
        phi_rim = S.phi(S.rim_idx);
        lambda_rim = S.lambda(S.rim_idx);
    case 'ideal'
        phi0  = 1.047198;
        theta = linspace(0, 2*pi, 400);
        phi_rim    = phi0 * ones(size(theta));   % latitude = phi0 everywhere
        lambda_rim = theta - pi;                 % or similar; longitude range choice
    otherwise
        error('Unknown rim_type "%s". Use "ideal" or "actual".', rim_type);
end
projection_reconstructedDataset(phi_data, lambda_data, ...
                                phi_rim, lambda_rim, ...
                                phi0, 'equidistant');

function [x, y, rho] = sphere_spherical_to_polar_cart(phi, lambda, preserve)
% sphere_spherical_to_polar_cart
% MATLAB analog of retistruct's sphere.spherical.to.polar.cart().
%
% Inputs
%   phi      : latitude (radians), column or row vector
%   lambda   : longitude (radians), same size as phi
%   preserve : 'latitude' (default) | 'area' | 'angle'
%              also accepts synonyms:
%                 equidistant -> latitude
%                 equalarea  -> area
%                 conformal / stereographic -> angle
%
% Outputs
%   x, y     : Cartesian coordinates on the polar disc (raw units)
%   rho      : radial coordinate prior to any rim scaling
%
% Notes
%   - Matches the south-polar aspect used in the Retistruct GUI.
%   - No scaling is applied to force the rim to unit radius; Retistruct
%     handles that later via rho.to.degrees()/phi.to.rho.

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

    % Retistruct uses x = rho*cos(lambda), y = rho*sin(lambda)
    x = rho .* cos(lambda);
    y = rho .* sin(lambda);
end


function projection_reconstructedDataset(phi_data, lambda_data, ...
                                         phi_rim,  lambda_rim, ...
                                         phi0, method)
% projection_reconstructedDataset
% MATLAB analog of retistruct::projection.reconstructedDataset()
%
% Inputs
%   phi_data, lambda_data : 1D vectors of point lat / lon (radians)
%   phi_rim,  lambda_rim  : outline rim lat / lon (radians)
%   phi0                  : rim latitude (same meaning as in Retistruct)
%   method                : 'equidistant' | 'equalarea' | 'conformal'
%
% Steps (matched to the R implementation):
%   1) project rim & datapoints to a polar disc (south-pole aspect)
%   2) convert the raw azimuthal coordinates into the GUI's degree grid
%   3) draw the circular retina, rim, and datapoints

    if nargin < 6
        method = 'equidistant';
    end

    [projectionMode, isAreaPreserving, projectionLabel] = normalize_projection(method);

    % Ensure column vectors
    phi_data    = phi_data(:);
    lambda_data = lambda_data(:);
    phi_rim     = phi_rim(:);
    lambda_rim  = lambda_rim(:);

    % --- 1. Project rim & datapoints to disc (raw units) ---
    [xr_raw, yr_raw] = sphere_spherical_to_polar_cart(phi_rim,  lambda_rim,  projectionMode);
    [xd_raw, yd_raw] = sphere_spherical_to_polar_cart(phi_data, lambda_data, projectionMode);

    % --- 2. Match Retistruct scaling (degrees with rim at phi0d+90) ---
    rim_xy   = rho_to_degrees([xr_raw, yr_raw], phi0, isAreaPreserving);
    data_xy  = rho_to_degrees([xd_raw, yd_raw], phi0, isAreaPreserving);
    rimRadius = phi_to_rho(phi0, phi0, isAreaPreserving);
    plotExtent = 1.05 * rimRadius;

    % Ideal rim circle for reference
    th  = linspace(0, 2*pi, 400);
    idealCircle = rimRadius * [cos(th); sin(th)];

    % --- 3. Plot setup ---
    figure('Color','w', 'Position', [100 100 1100 500]);

    % Subplot 1: scatter + rim (Retistruct style)
    subplot(1,2,1); hold on; axis equal;
    set(gca, 'Visible', 'off');
    plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1);
    plot(rim_xy(:,1), rim_xy(:,2), 'k-', 'LineWidth', 1.5);
    plot(data_xy(:,1), data_xy(:,2), 'o', ...
        'MarkerSize', 2, ...
        'MarkerFaceColor', [0.2 0.4 1.0], ...
        'MarkerEdgeColor', 'k');
    xlim([-plotExtent, plotExtent]);
    ylim([-plotExtent, plotExtent]);
    title(sprintf('Retistruct-style projection (%s)', projectionLabel), 'Interpreter', 'none');

    % Subplot 2: adaptive KDE density heatmap on flattened disc
    subplot(1,2,2); hold on; axis equal;
    gridRes = 220;
    gridLin = linspace(-plotExtent, plotExtent, gridRes);
    [gridX, gridY] = meshgrid(gridLin, gridLin);
    insideMask = hypot(gridX, gridY) <= (rimRadius + 1e-6);
    densityGrid = nan(size(gridX));
    if any(insideMask(:))
       evalPoints = [gridX(insideMask), gridY(insideMask)];
       densityVals = adaptiveKDE(data_xy, evalPoints);
       densityGrid(insideMask) = densityVals;
    end
    contourf(gridX, gridY, densityGrid, 32, 'LineColor', 'none');
    colormap(gca, parula);
    plot(idealCircle(1,:), idealCircle(2,:), 'k-', 'LineWidth', 1.25);
    plot(rim_xy(:,1), rim_xy(:,2), 'k-', 'LineWidth', 1.5);
    xlim([-plotExtent, plotExtent]);
    ylim([-plotExtent, plotExtent]);
    axis off;
    title('Adaptive KDE density on flattened disc');
    cb = colorbar('Location', 'eastoutside');
    cb.Label.String = 'density';
end

function densityVals = adaptiveKDE(points, evalPoints)
% adaptiveKDE - Adaptive Gaussian KDE with Abramson's method
%
% points:     [n x 2] sample points (x,y)
% evalPoints: [m x 2] evaluation points (x,y)
% densityVals: [m x 1] density estimates

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
% Scale polar Cartesian coordinates to match Retistruct's degree grid.
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
% Convert a latitude to the GUI's radial coordinate.
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
% Area-preserving azimuthal radius used by Retistruct (unit sphere).
    rho = sqrt(2 * (1 + sin(phi)));
end

function [mode, isAreaPreserving, label] = normalize_projection(method)
% Map user-friendly method names to Retistruct projection modes.
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

