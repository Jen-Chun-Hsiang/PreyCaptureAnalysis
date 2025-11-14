%% 1) CIRCULAR VISUAL FIELD DISK - Azimuthal Projection centered at (0°, 0°)
% For publication-quality visual field maps (like panels A₃/C₃):
%   - Use Lambert azimuthal equal-area projection (area-preserving, best for density)
%   - OR stereographic projection (conformal, preserves local shapes)
%   - Centered at visual field origin (0° azimuth, 0° elevation)
%
% YOUR DATA FORMAT:
%   x = elevation (lambda) in radians: -π/2 to +π/2 (-90° to +90°)
%   y = azimuth (phi) in radians: -π to +π (-180° to +180°)

dataFile = 'C:\Users\jhsiang\Documents\R-retistruct\Test06\data.mat';



% DIAGNOSTIC: Inspect data ranges and distribution
S = load(dataFile);
fprintf('=== DATA INSPECTION ===\n');
fprintf('x (elevation) range: [%.4f, %.4f] rad = [%.1f°, %.1f°]\n', ...
    min(S.x), max(S.x), rad2deg(min(S.x)), rad2deg(max(S.x)));
fprintf('y (azimuth) range: [%.4f, %.4f] rad = [%.1f°, %.1f°]\n', ...
    min(S.y), max(S.y), rad2deg(min(S.y)), rad2deg(max(S.y)));
fprintf('Number of points: %d\n', length(S.x));

%%
clc; close all;
x = S.x; y = S.y/2;
shift_x = (max(x) + min(x))/2;
x = (x-shift_x)*pi/range(x);
X = x(:);Y = y(:);

% x = (-90:10:0).*(pi/180);   % fixed azimuth
% y = (-90:10:90).*(pi/180);               % elevations
% [X, Y] = meshgrid(x, y);
% X = X(:);Y = Y(:);
figure; subplot(2,2,1);
scatter(X, Y, 15, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('y (azimuth)'); ylabel('x (elevation)');
subplot(2,2,2);
rho = acos( cos(Y) .* cos(X) );
Xt = cos(Y) .* sin(X);
Yt = sin(Y);
theta = atan2(Yt, Xt);
% polarplot(theta, rho);
polarscatter(theta, rho,  15 , 'filled');

nbins = 15;

xedges = linspace(min(X), max(X), nbins+1);
yedges = linspace(min(Y), max(Y), nbins+1);

[N, xedges, yedges] = histcounts2(X, Y, xedges, yedges);

% Bin centers
xc = (xedges(1:end-1) + xedges(2:end)) / 2;
yc = (yedges(1:end-1) + yedges(2:end)) / 2;
[AZgrid, ELgrid] = meshgrid(xc, yc);  % same coordinates as density N

subplot(2,2,3);
imagesc(xc, yc, N');      % note transpose so axes match
set(gca, 'YDir', 'normal');
hold on;
contour(AZgrid, ELgrid, N', 6, 'k', 'LineWidth', 1); % 6 contour levels
xlabel('azimuth (rad)');
ylabel('elevation (rad)');
title('Density contour in (az, el)');
axis tight; axis equal; grid on;
colorbar;

% ---- Point projection (your original logic) ----
rho_pts = acos( cos(Y) .* cos(X) );
Xp = cos(Y) .* sin(X);
Yp = sin(Y);
theta_pts = atan2(Yp, Xp);

% ---- Grid projection for the contour ----
rho_grid = acos( cos(ELgrid) .* cos(AZgrid) );
Xg = cos(ELgrid) .* sin(AZgrid);
Yg = sin(ELgrid);
theta_grid = atan2(Yg, Xg);

% Convert polar (theta, rho) grid to Cartesian for plotting
Xproj = rho_grid .* cos(theta_grid);
Yproj = rho_grid .* sin(theta_grid);

% Convert the point cloud too (for overlay)
Xpts_proj = rho_pts .* cos(theta_pts);
Ypts_proj = rho_pts .* sin(theta_pts);

subplot(2,2,4);
hold on;

% Projected contour (same density N, now warped by projection)
contour(Xproj, Yproj, N', 6, 'k', 'LineWidth', 1);  % same 6 levels

% Projected scatter
scatter(Xpts_proj, Ypts_proj, 15, 'filled', 'MarkerFaceAlpha', 0.4);

% Make it look like a polar map
axis equal;
xlabel('x (polar)');
ylabel('y (polar)');
title('Projected contour in polar map');
grid on;

% Optional: draw outer circle if you like
maxR = max(rho_pts(:));
t = linspace(0, 2*pi, 400);
plot(maxR*cos(t), maxR*sin(t), 'r--');

%%
% X, Y are your azimuth/elevation in radians
az = X(:);
el = Y(:);
points = [az, el];

% Grid for KDE evaluation
ngrid = 100;  % increase for smoother, decrease for speed
az_lin = linspace(min(az), max(az), ngrid);
el_lin = linspace(min(el), max(el), ngrid);
[AZgrid, ELgrid] = meshgrid(az_lin, el_lin);
gridPoints = [AZgrid(:), ELgrid(:)];

% Adaptive KDE on (az, el)
densityVals = adaptiveKDE(points, gridPoints);
D = reshape(densityVals, size(AZgrid));   % D(i,j) corresponds to (AZgrid(i,j), ELgrid(i,j))

figure;

% 1) Scatter in (az, el) space
subplot(1,3,1);
scatter(az, el, 15, 'filled', 'MarkerFaceAlpha', 0.4);
xlabel('azimuth (rad)');
ylabel('elevation (rad)');
title('Raw samples (az, el)');
axis equal; grid on;

% 2) Smooth KDE contour in (az, el) space
subplot(1,3,2);
contourf(AZgrid, ELgrid, D, 12, 'LineColor', 'none');  % filled smooth map
hold on;
contour(AZgrid, ELgrid, D, 8, 'k', 'LineWidth', 0.7);  % contour lines
scatter(az, el, 8, 'k.', 'MarkerFaceAlpha', 0.2);      % optional: show points
xlabel('azimuth (rad)');
ylabel('elevation (rad)');
title('Adaptive KDE in (az, el)');
axis equal; grid on;
colorbar;

% 3) Project everything to polar (theta, rho)

% --- Project points ---
rho_pts = acos( cos(el) .* cos(az) );
Xp = cos(el) .* sin(az);
Yp = sin(el);
theta_pts = atan2(Yp, Xp);

% --- Project grid (same formula applied to grid) ---
rho_grid = acos( cos(ELgrid) .* cos(AZgrid) );
Xg = cos(ELgrid) .* sin(AZgrid);
Yg = sin(ELgrid);
theta_grid = atan2(Yg, Xg);

% Convert to Cartesian for plotting (so we can contour on a regular axes)
Xproj = rho_grid .* cos(theta_grid);
Yproj = rho_grid .* sin(theta_grid);

% And the projected points
Xpts_proj = rho_pts .* cos(theta_pts);
Ypts_proj = rho_pts .* sin(theta_pts);

subplot(1,3,3);
hold on;

% Smooth contour map in polar projection
contourf(Xproj, Yproj, D, 12, 'LineColor', 'none');
contour(Xproj, Yproj, D, 8, 'k', 'LineWidth', 0.7);

% Projected scatter on top
scatter(Xpts_proj, Ypts_proj, 10, 'k.', 'MarkerFaceAlpha', 0.3);

axis equal;
xlabel('x (polar)');
ylabel('y (polar)');
title('Adaptive KDE contour in polar projection');
grid on;
colorbar;

%%
keyboard
%%

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



%%
% Calculate optimal R based on data extent
% First project with a provisional R to determine actual extent
opts_temp = struct('projection','lambert', 'angleUnit','rad', ...
                   'centerPhi',0, 'centerLambda',0, ...
                   'R',sqrt(2), ...
                   'drawGrid',true, 'gridStepDeg',30, ...
                   'markerSize', 8);
[X_temp, Y_temp] = project_visual_field(S.x, S.y, opts_temp);
max_radius = max(hypot(X_temp, Y_temp));

% Now set R based on actual data extent
opts = struct('projection','lambert', 'angleUnit','rad', ...
              'centerPhi',0, 'centerLambda',0, ...
              'R',max_radius, ...
              'drawGrid',true, 'gridStepDeg',30, ...
              'markerSize', 8);

[X2D, Y2D] = project_visual_field(S.x, S.y, opts);
figure; hold on; axis equal off
th = linspace(0,2*pi,720); 
plot(max_radius*cos(th), max_radius*sin(th), 'k-'); % disk boundary based on data extent
scatter(X2D, Y2D, 12, 'filled');

% Check if values look like they could be Cartesian (not angles)
if max(abs(S.x)) < 0.1 && max(abs(S.y)) < 0.1
    warning('Values are very small (< 0.1). Are these mm coordinates, not angles?');
elseif max(abs(S.x)) > pi/2 || max(abs(S.y)) > pi
    warning('Values exceed expected angular ranges. Check coordinate system.');
end

keyboard
% DIAGNOSTIC PLOTS: Compare different interpretations
figure('Color','w', 'Position', [100 100 1400 500]);

% Plot 1: Direct scatter (treating as Cartesian-like)
subplot(1,3,1);
scatter(S.y, S.x, 15, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('y (azimuth)'); ylabel('x (elevation)');
title('Raw Data (as-is)');
axis equal; grid on;

% Plot 2: Converted to degrees (assuming radians)
subplot(1,3,2);
scatter(rad2deg(S.y), rad2deg(S.x), 15, 'filled', 'MarkerFaceAlpha', 0.5);
xlabel('Azimuth (degrees)'); ylabel('Elevation (degrees)');
title('Interpreted as Radians → Degrees');
axis equal; grid on;

% Plot 3: Azimuthal projection
subplot(1,3,3);
opts = struct('projection','lambert', 'angleUnit','rad', ...
              'centerPhi',0, 'centerLambda',0, ...
              'R',sqrt(2), ...
              'drawGrid',true, 'gridStepDeg',30, ...
              'markerSize', 8);
[X2D, Y2D] = plot_visual_field_from_mat_simple(dataFile, opts);
title('Lambert Projection');

% Print diagnostic info about projection
fprintf('\n=== PROJECTION DIAGNOSTIC ===\n');
fprintf('Projected radial distances: min=%.4f, max=%.4f, mean=%.4f\n', ...
    min(hypot(X2D, Y2D)), max(hypot(X2D, Y2D)), mean(hypot(X2D, Y2D)));
fprintf('Points near edge (r > 0.9*R): %d / %d (%.1f%%)\n', ...
    sum(hypot(X2D, Y2D) > 0.9*sqrt(2)), length(X2D), ...
    100*sum(hypot(X2D, Y2D) > 0.9*sqrt(2))/length(X2D));

%%
% 2) If you already have phi/lambda in workspace (vectors):
phi = randn(100,1)*30;      % deg, example azimuths
lambda = randn(100,1)*20;   % deg, example elevations
[X2D, Y2D] = project_visual_field(phi, lambda, opts);  % and then plot as you wish
figure; scatter(X2D, Y2D, 18, 'filled'); axis equal off;

%%
% Synthetic phi (azimuth) and lambda (elevation) data
% Arrange points in a spherical grid: -60° to +60° azimuth, -45° to +45° elevation
[phi_grid, lambda_grid] = meshgrid(-60:15:60, -45:15:45);

% Flatten to vectors
phi = phi_grid(:);
lambda = lambda_grid(:);

% Optional: save to a .mat file to mimic your R-exported data
x = phi;
y = lambda;
save('synthetic_phi_lambda.mat', 'x', 'y');

% Set up plotting options
opts = struct( ...
    'projection', 'lambert', ...
    'angleUnit', 'deg', ...
    'centerPhi', 0, ...
    'centerLambda', 0, ...
    'R', 1, ...
    'drawGrid', true, ...
    'gridStepDeg', 30 ...
);

% Call the plotting function
[X2D, Y2D] = plot_visual_field_from_mat('synthetic_phi_lambda.mat', opts);

%%

function [X2D, Y2D] = plot_visual_field_from_mat(matFile, opts)
% plot_visual_field_from_mat Load elevation/azimuth and create visual field map
% For RETINOTOPIC MAPPING: projects visual angles onto 2D disk using spherical projection
%
% Inputs:
%   matFile : string path to .mat file (with variables 'x' and 'y')
%             x = elevation (lambda) in radians or degrees
%             y = azimuth (phi) in radians or degrees
%   opts    : struct with fields (all optional):
%               .projection   : 'lambert' (equal-area, default) or 'stereo'
%               .angleUnit    : 'rad' (default) or 'deg'
%               .centerPhi    : center azimuth (default 0)
%               .centerLambda : center elevation (default 0)
%               .R            : disk radius (default 1)
%               .drawGrid     : true/false (default true)
%               .gridStepDeg  : spacing for grid (default 30)
%               .markerSize   : default 18
%               .markerColor  : default [0 0.45 0.74]
%
% Outputs:
%   X2D, Y2D : projected 2D coordinates

    S = load(matFile);
    if ~isfield(S, 'x') || ~isfield(S, 'y')
        error('Expected variables x and y in %s', matFile);
    end
    % YOUR DATA FORMAT: x = elevation, y = azimuth (both in radians)
    lambda = S.x;    % elevation (latitude) in radians
    phi = S.y;       % azimuth (longitude) in radians

    if nargin < 2 || isempty(opts), opts = struct(); end
    if ~isfield(opts, 'projection'),   opts.projection   = 'lambert'; end
    if ~isfield(opts, 'angleUnit'),    opts.angleUnit    = 'deg'; end
    if ~isfield(opts, 'centerPhi'),    opts.centerPhi    = 0; end
    if ~isfield(opts, 'centerLambda'), opts.centerLambda = 0; end
    if ~isfield(opts, 'R'),            opts.R            = 1; end
    if ~isfield(opts, 'drawGrid'),     opts.drawGrid     = true; end
    if ~isfield(opts, 'gridStepDeg'),  opts.gridStepDeg  = 30; end
    if ~isfield(opts, 'markerSize'),   opts.markerSize   = 18; end
    if ~isfield(opts, 'markerColor'),  opts.markerColor  = [0 0.45 0.74]; end

    [X2D, Y2D] = project_visual_field(phi, lambda, opts);

    % ---- Plot ----
    figure('Color','w'); hold on;
    % Disk boundary
    th = linspace(0, 2*pi, 720);
    plot(opts.R*cos(th), opts.R*sin(th), 'k-', 'LineWidth', 1.25);

    % Optional grid (meridians/parallels)
    if opts.drawGrid
        draw_visual_field_grid(opts);
    end

    % Points
    scatter(X2D, Y2D, opts.markerSize, 'filled', 'MarkerFaceColor', opts.markerColor, ...
            'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', 'none');

    axis equal; xlim(opts.R*[-1 1]); ylim(opts.R*[-1 1]);
    set(gca, 'Visible', 'off');
    title(sprintf('Visual Field (%s projection, center=(%g,%g) %s)', ...
          opts.projection, opts.centerPhi, opts.centerLambda, opts.angleUnit), ...
          'FontWeight','normal');
end


function [X2D, Y2D] = project_visual_field(phi, lambda, opts)
% project_visual_field Project spherical points (phi, lambda) to a 2D disk.
% phi      : azimuth (longitude)
% lambda   : elevation (latitude)
% opts     : same struct as above (angleUnit, projection, center, R)

    % Units → radians
    if strcmpi(opts.angleUnit, 'deg')
        d2r = pi/180;
        phi = phi * d2r;
        lambda = lambda * d2r;
        phi0 = opts.centerPhi * d2r;
        lam0 = opts.centerLambda * d2r;
    else
        phi0 = opts.centerPhi;
        lam0 = opts.centerLambda;
    end

    % Δφ relative to center (wrap to [-pi, pi] for numerical stability)
    dphi = wrapToPi(phi - phi0);

    % Great-circle angular distance from center:
    % cosθ = sin(lam0) sin(λ) + cos(lam0) cos(λ) cos(Δφ)
    cosTheta = sin(lam0).*sin(lambda) + cos(lam0).*cos(lambda).*cos(dphi);
    cosTheta = min(max(cosTheta, -1), 1); % clamp
    theta = acos(cosTheta);

    % Bearing α from center to point (for disk angle)
    % α = atan2( cos(λ)*sin(Δφ), cos(lam0)*sin(λ) - sin(lam0)*cos(λ)*cos(Δφ) )
    alpha = atan2( cos(lambda).*sin(dphi), ...
                   cos(lam0).*sin(lambda) - sin(lam0).*cos(lambda).*cos(dphi) );

    % Radial mapping r(θ)
    R = opts.R;
    switch lower(opts.projection)
        case 'lambert'     % Lambert azimuthal equal-area
            r = R .* sqrt(2) .* sin(theta ./ 2);  % r = R * sqrt(2) * sin(θ/2)
        case 'stereo'      % Stereographic
            r = 2 .* R .* tan(theta ./ 2);        % r = 2R * tan(θ/2)
        otherwise
            error('Unknown projection "%s". Use ''lambert'' or ''stereo''.', opts.projection);
    end

    X2D = r .* cos(alpha);
    Y2D = r .* sin(alpha);

    % Clip tiny FP noise just outside disk
    clip = hypot(X2D, Y2D) > (R + 1e-10);
    X2D(clip) = X2D(clip) .* (R ./ hypot(X2D(clip), Y2D(clip)));
    Y2D(clip) = Y2D(clip) .* (R ./ hypot(X2D(clip), Y2D(clip)));
end


function draw_visual_field_grid(opts)
% draw_visual_field_grid Draw meridians/parallels onto the disk (visual aid).
    R = opts.R;
    step = opts.gridStepDeg;
    th = linspace(0, 2*pi, 720);

    % Light boundary (already drawn bold once)
    plot(R*cos(th), R*sin(th), 'Color', 0.8*[1 1 1]);

    % Concentric rings for zenith angle (every step degrees from center)
    hold on;
    for ang = step:step:(180 - step)
        % For Lambert, r = R * sqrt(2) * sin(θ/2); For Stereo, r = 2R * tan(θ/2)
        theta = deg2rad(ang);
        switch lower(opts.projection)
            case 'lambert'
                r = R * sqrt(2) * sin(theta/2);
            case 'stereo'
                r = 2*R * tan(theta/2);
                r = min(r, R); % cap to disk
        end
        if r < R - 1e-6
            plot(r*cos(th), r*sin(th), '-', 'Color', 0.9*[1 1 1]);
        end
    end

    % Spokes (meridians) every step degrees
    for az = 0:step:(360-step)
        a = deg2rad(az);
        plot([0, R*cos(a)], [0, R*sin(a)], '-', 'Color', 0.9*[1 1 1]);
    end
end


function ang = wrapToPi(ang)
% wrapToPi Wrap angle in radians to [-pi, pi]
    ang = mod(ang + pi, 2*pi) - pi;
end


function [X2D, Y2D] = plot_visual_field_from_mat_simple(matFile, opts)
% Simplified version without creating new figure - for subplots
    S = load(matFile);
    lambda = S.x;    % elevation
    phi = S.y;       % azimuth
    
    if ~isfield(opts, 'markerSize'), opts.markerSize = 18; end
    if ~isfield(opts, 'markerColor'), opts.markerColor = [0 0.45 0.74]; end
    
    [X2D, Y2D] = project_visual_field(phi, lambda, opts);
    
    % Plot on current axes
    hold on;
    th = linspace(0, 2*pi, 720);
    plot(opts.R*cos(th), opts.R*sin(th), 'k-', 'LineWidth', 1.25);
    
    if opts.drawGrid
        draw_visual_field_grid(opts);
    end
    
    scatter(X2D, Y2D, opts.markerSize, 'filled', 'MarkerFaceColor', opts.markerColor, ...
            'MarkerFaceAlpha', 0.5, 'MarkerEdgeColor', 'none');
    
    axis equal; xlim(opts.R*[-1 1]); ylim(opts.R*[-1 1]);
    set(gca, 'Visible', 'off');
end
