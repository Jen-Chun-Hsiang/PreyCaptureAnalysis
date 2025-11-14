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

function [x, y, rho] = sphere_spherical_to_polar_cart(phi, lambda, phi0, method)
% sphere_spherical_to_polar_cart
% MATLAB analog of retistruct's sphere.spherical.to.polar.cart()
% plus azimuthal.*() projections (south-polar aspect).
%
% Inputs
%   phi    : latitude (radians), column or row vector
%   lambda : longitude (radians), same size as phi
%   phi0   : rim latitude (radians) – the retinal edge latitude
%   method : 'equidistant' | 'equalarea' | 'conformal'
%            (matches azimuthal.equidistant / equalarea / conformal)
%
% Outputs
%   x, y   : Cartesian coordinates on the polar disc
%   rho    : radial coordinate (before scaling to unit rim)
%
% All formulas use a unit sphere and *south-polar* aspect:
%   - South pole → centre of disc
%   - Rim at latitude phi0 → unit circle

    if nargin < 4
        method = 'equidistant';  % Retistruct often uses equidistant
    end

    % Ensure column vectors
    phi    = phi(:);
    lambda = lambda(:);

    % --- 1. Convert latitude to colatitude ψ from NORTH pole ---
    % Standard spherical: ψ = 0 at north pole, π at south pole
    psi = pi/2 - phi;        % ψ in [0, π]

    % --- 2. Raw radial coordinate ρ_raw depending on projection ---
    method = lower(method);
    switch method
        case 'equidistant'
            % Azimuthal equidistant, centred on SOUTH pole
            % Central angle from south pole: c = π - ψ
            % ρ ∝ c
            rho_raw = pi - psi;  % = π/2 + phi

        case 'equalarea'
            % Lambert azimuthal equal-area, south-polar aspect
            % From Wikipedia / MathWorld in colatitude ψ (north-based):
            %   R = 2 * cos(ψ/2) for centre at south pole. :contentReference[oaicite:1]{index=1}
            rho_raw = 2 * cos(psi / 2);

        case 'conformal'
            % Stereographic / azimuthal conformal, south-polar aspect
            % Central angle from south pole: c = π - ψ
            % Stereographic radius for central angle c:
            %   R = 2 * tan(c/2). :contentReference[oaicite:2]{index=2}
            c       = pi - psi;
            rho_raw = 2 * tan(c / 2);

        otherwise
            error('Unknown method "%s". Use equidistant | equalarea | conformal.', method);
    end

    % --- 3. Scale so that rim latitude phi0 → unit circle (ρ = 1) ---
    psi0 = pi/2 - phi0;   % colatitude of rim

    switch method
        case 'equidistant'
            rho0 = pi - psi0;              % = π/2 + phi0
        case 'equalarea'
            rho0 = 2 * cos(psi0 / 2);
        case 'conformal'
            c0   = pi - psi0;
            rho0 = 2 * tan(c0 / 2);
    end

    scale = 1 / rho0;
    rho   = scale * rho_raw;

    % --- 4. Convert to Cartesian in polar plane ---
    % Choice of angular convention: here we take
    %   λ = 0 → +Y axis
    %   λ increases counter-clockwise
    % This is:
    %   x = ρ * sin(λ)
    %   y = ρ * cos(λ)
    %
    % If you find the plot rotated vs. the Retistruct GUI,
    % you can swap sin/cos or flip signs later.
    x = rho .* sin(lambda);
    y = rho .* cos(lambda);

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
% This:
%   1) projects rim & datapoints to a polar disc
%   2) draws the circular retina
%   3) overlays the datapoints

    if nargin < 6
        method = 'equidistant';
    end

    % Ensure column vectors
    phi_data    = phi_data(:);
    lambda_data = lambda_data(:);
    phi_rim     = phi_rim(:);
    lambda_rim  = lambda_rim(:);

    % --- 1. Project rim & datapoints to disc ---
    [xr, yr] = sphere_spherical_to_polar_cart(phi_rim,  lambda_rim,  phi0, method);
    [xd, yd] = sphere_spherical_to_polar_cart(phi_data, lambda_data, phi0, method);

    % --- 2. Plot setup ---
    figure('Color','w', 'Position', [100 100 1100 500]);

    % Subplot 1: original Retistruct-style scatter + rim
    subplot(1,2,1); hold on; axis equal;
    set(gca, 'Visible', 'off');  % like Retistruct’s default

    % Draw rim as a smooth circle at radius 1 (idealised boundary)
    th  = linspace(0, 2*pi, 400);
    xc  = cos(th);
    yc  = sin(th);
    plot(xc, yc, 'k-', 'LineWidth', 1);    % ideal rim
    % Also draw the actual rim polygon from the reconstruction
    plot(xr, yr, 'k-', 'LineWidth', 1.5);
    % --- 3. Plot datapoints ---
    plot(xd, yd, 'o', ...
        'MarkerSize', 2, ...
        'MarkerFaceColor', [0.2 0.4 1.0], ...
        'MarkerEdgeColor', 'k');
    xlim([-1.1, 1.1]);
    ylim([-1.1, 1.1]);
    title(sprintf('Retistruct-style projection (%s)', method), 'Interpreter', 'none');

    % Subplot 2: adaptive KDE density heatmap on flattened disc
    subplot(1,2,2); hold on; axis equal;
    gridLimit = 1.05;
    gridRes = 220;
    gridLin = linspace(-gridLimit, gridLimit, gridRes);
    [gridX, gridY] = meshgrid(gridLin, gridLin);
    insideMask = hypot(gridX, gridY) <= 1 + 1e-6;
    densityGrid = nan(size(gridX));
    if any(insideMask(:))
       evalPoints = [gridX(insideMask), gridY(insideMask)];
       densityVals = adaptiveKDE([xd, yd], evalPoints);
       densityGrid(insideMask) = densityVals;
    end
    contourf(gridX, gridY, densityGrid, 32, 'LineColor', 'none');
    colormap(gca, parula);
    plot(xc, yc, 'k-', 'LineWidth', 1.25);
    plot(xr, yr, 'k-', 'LineWidth', 1.5);
    xlim([-gridLimit, gridLimit]);
    ylim([-gridLimit, gridLimit]);
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

