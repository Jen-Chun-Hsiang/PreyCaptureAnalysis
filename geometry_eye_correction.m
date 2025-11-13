% Binocular RF tilt correction on a spherical dome
% Avah's setup:
%   - Spherical dome radius (center-to-canvas): R = 24 cm
%   - Inter-ocular distance: d = 1 cm
%   - Tilt gamma: -5:1:5 degrees (right eye closer when gamma > 0)
%   - RF sets (measured canvas azimuths, degrees):
%       Set 1: L = -3, R = +3
%       Set 2: L = -3, R = +2
%       Set 3: L = -3, R = +1
%
% Output:
%   - Plots the binocular RF shift (true - observed) vs tilt angle
%     for the three sets.

clear; clc;

% Geometry
R = 24;       % cm, radius of the dome (centered at midpoint between eyes)
d = 1;        % cm, inter-ocular distance
gammas = -5:1:5;   % deg, right eye is closer for positive values

% Measured (canvas) azimuths in degrees: [theta_L, theta_R]
sets = [
    -3,  3;   % Set 1
    -3,  2;   % Set 2
    -3,  1    % Set 3
];

% Helper: wrap to [-180, 180)
wrap180 = @(a) mod(a + 180, 360) - 180;

% Given a canvas azimuth theta (deg), return the 3D point on the sphere
spherePoint = @(theta_deg) [ R * sind(theta_deg), 0, R * cosd(theta_deg) ];

% Compute the "true" azimuth (deg) seen from an eye at position E (cm)
% relative to the head-centered z-axis. We take the azimuth of vector S - E.
trueAzFromEye = @(theta_obs_deg, E) atan2d( ...
    (R * sind(theta_obs_deg) - E(1)), ...
    (R * cosd(theta_obs_deg) - E(3)) );

% Preallocate: each row corresponds to a set; columns across gamma
nSets   = size(sets,1);
nGamma  = numel(gammas);
shift_deg = zeros(nSets, nGamma);   % (true separation) - (observed separation)
true_sep  = zeros(nSets, nGamma);   % binocular separation after correction
obs_sep   = zeros(nSets, 1);        % constant per set

for s = 1:nSets
    thetaL = sets(s,1);  % measured canvas azimuth (deg) for left eye
    thetaR = sets(s,2);  % measured canvas azimuth (deg) for right eye
    obs_sep(s) = wrap180(thetaR - thetaL);

    for i = 1:nGamma
        g = gammas(i);

        % Depth asymmetry: right eye closer for g > 0 (front-back axis = +z)
        % We model only the front-back asymmetry (small angles), keeping lateral ±d/2.
        % Depth offset between eyes ~ d * sind(g); split half to each eye.
        dz = (d/2) * sind(g);
        ER = [ +d/2, 0, +dz];   % Right eye position (cm)
        EL = [ -d/2, 0, -dz];   % Left  eye position (cm)

        % True azimuths seen from each eye to the same canvas points
        % Note: inline evaluate spherePoint; avoid function-handle indexing issues
        SR = [R*sind(thetaR), 0, R*cosd(thetaR)];
        SL = [R*sind(thetaL), 0, R*cosd(thetaL)];

        psiR = atan2d(SR(1) - ER(1), SR(3) - ER(3));  % true azimuth from right eye
        psiL = atan2d(SL(1) - EL(1), SL(3) - EL(3));  % true azimuth from left eye

        % Binocular separation (true) and induced shift
        true_sep(s,i) = wrap180(psiR - psiL);
        shift_deg(s,i) = wrap180(true_sep(s,i) - obs_sep(s));
    end
end

% ---- Plot: induced shift vs tilt angle for each set ----
figure('Color','w'); hold on;
plot(gammas, shift_deg(1,:), 'LineWidth', 2);
plot(gammas, shift_deg(2,:), 'LineWidth', 2);
plot(gammas, shift_deg(3,:), 'LineWidth', 2);
yline(0,':','Color',[0.4 0.4 0.4]);

xlabel('Tilt angle \gamma (deg)  [right eye closer if \gamma > 0]');
ylabel('RF shift (true - observed) (deg)');
title('Induced binocular RF shift vs inter-eye tilt on a 24 cm dome (d = 1 cm)');
legend({'Set 1: (-3, +3)', 'Set 2: (-3, +2)', 'Set 3: (-3, +1)'}, 'Location','best');
grid on;

% (Optional) Also plot the corrected binocular separations
figure('Color','w'); hold on;
plot(gammas, true_sep(1,:), 'LineWidth', 2);
plot(gammas, true_sep(2,:), 'LineWidth', 2);
plot(gammas, true_sep(3,:), 'LineWidth', 2);
yline(obs_sep(1), ':', 'Color', [0.6 0.6 0.6], 'Label', 'Set 1 observed');
yline(obs_sep(2), ':', 'Color', [0.6 0.6 0.6], 'Label', 'Set 2 observed');
yline(obs_sep(3), ':', 'Color', [0.6 0.6 0.6], 'Label', 'Set 3 observed');

xlabel('Tilt angle \gamma (deg)  [right eye closer if \gamma > 0]');
ylabel('Binocular separation (deg)');
title('Corrected (true) binocular separations vs tilt');
legend({'Set 1 true','Set 2 true','Set 3 true'}, 'Location','best');
grid on;

% Print a quick table at a few tilt angles for sanity
probeGammas = [-5 0 5];
fprintf('\nSummary (true separation in deg) at gamma = -5, 0, +5:\n');
for s = 1:nSets
    [~,idx] = ismember(probeGammas, gammas);
    vals = true_sep(s, idx);
    fprintf('  Set %d  (obs = %.1f):   %.3f   %.3f   %.3f\n', ...
        s, obs_sep(s), vals(1), vals(2), vals(3));
end
