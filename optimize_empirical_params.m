function [fac_opt, offset_opt] = optimize_empirical_params(rgc_time_std, equivalent_noise_level, nRepeats, r_target, mean_rt, rect_thr, T, initial_guess)
%OPTIMIZE_EMPIRICAL_PARAMS  Fit noise scale (fac) and baseline offset
%   such that simulated reliability matches r_target and the mean post-
%   rectification equals mean_rt.
%
%  [fac_opt, offset_opt] = optimize_empirical_params(
%       rgc_time_std,           % STD of smoothed, rectified rgc_time fluctuations
%       equivalent_noise_level, % pre-smoothing noise STD x
%       nRepeats,               % number of repeats to simulate
%       r_target,               % desired repeat-reliability
%       mean_rt,                % observed mean of smoothed, rectified rgc_time
%       rect_thr,               % rectification threshold used in sim
%       T,                      % number of timepoints per trace
%       initial_guess           % [fac0, offset0] initial guesses
%  )

    % Set fsolve options
    opts = optimset('Display','off', 'TolFun',1e-6, 'TolX',1e-6);
    % Solve for both parameters
    sol = fsolve(@(x) error_funcs(x, rgc_time_std, equivalent_noise_level, nRepeats, r_target, mean_rt, rect_thr, T), initial_guess, opts);

    fac_opt    = sol(1);
    offset_opt = sol(2);
end

function F = error_funcs(x, rgc_time_std, equivalent_noise_level, nRepeats, r_target, mean_rt, rect_thr, T)
%ERROR_FUNCS  Compute errors for reliability and mean constraints
%   x(1) = fac, x(2) = offset

    fac    = x(1);
    offset = x(2);
    nTr = nRepeats * 10;

    % 1) Generate raw traces around offset
    raw = offset + rgc_time_std * randn(nTr, T);

    % 2) Add pre-rectification (white) noise scaled by fac
    noisy_raw = raw + equivalent_noise_level * fac * randn(size(raw));

    % 3) Apply rectification at rect_thr
    rectified = max(noisy_raw, rect_thr);

    % 4) Compute mean after rectification
    mean_sim = mean(rectified(:));

    % 5) Compute reliability r_s from pairwise correlations
    R = corr(rectified');
    mask = triu(true(size(R)), 1);
    r_s  = mean(R(mask));

    % 6) Errors: [r_s - r_target; mean_sim - mean_rt]
    F = [r_s - r_target;
         mean_sim - mean_rt];
end
