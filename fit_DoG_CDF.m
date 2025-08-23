
% Results struct is left in workspace as "results"

% ---------------- Enhanced DoG-CDF Fitting Function ----------------
function [s_plot, y_plot, p_struct, stats] = fit_DoG_CDF(s, y)
% Fits Difference-of-Gaussians CDF model (center-surround) using standard CDFs:
% y(s) = global_bias + gain * (CDF_center(s) + w_s * CDF_surround(s))
%
% Notes for initialization (matching your example intuition):
% - x starts near the center of the CDF (e.g., x=0, mu_c≈0 => CDF≈0.5)
% - Surround has slightly larger mean (mu_s > mu_c) and wider sigma (≈4x)
% - Surround weight starts small and subtractive (w_s ≈ -0.1)
%
% Constraints:
%   - Center: mu_c, sigma_c > 0
%   - Surround: mu_s >= mu_c, sigma_s >= sigma_c (broader/larger)
%   - gain >= 0 (positive scaling)
%   - w_s in (-1, 0] (subtractive surround)
%   - global_bias is free

s = s(:); y = y(:);
if length(s) ~= length(y) || length(s) < 3
    error('Invalid input: s and y must be same length with at least 3 points');
end

% Data normalization for numerical stability
s_std = std(s);
if s_std == 0, s_std = 1; end
s_norm = s / s_std;

y_mean = mean(y); y_std = std(y);
if y_std == 0, y_std = 1; end
y_norm = (y - y_mean) / y_std;

% Initial parameter guesses (in normalized space)
s_min = min(s_norm); s_max = max(s_norm); s_range = s_max - s_min;
y_range = max(y_norm) - min(y_norm);

% Initialize 6 parameters: [mu_c, log_sigma_c, delta_mu_s, log_delta_sigma_s, log_gain, w_s]
p0 = zeros(6,1);
% Center mean near the smallest size (so starting around the CDF midpoint)
p0(1) = s_min + 0.05 * s_range;      % mu_c (center position)
% Center width as ~20% of range
p0(2) = log( max(0.2 * s_range, 1e-3) );         % log(sigma_c)
% Surround mean a bit larger than center
p0(3) = 0.3 * s_range;               % delta_mu_s = mu_s - mu_c >= 0
% Surround width about 4x center: delta_sigma ≈ 3*sigma_c
sigma_c_init = exp(p0(2));
p0(4) = log( max(3 * sigma_c_init, 1e-3) );      % log(delta_sigma_s) = log(sigma_s - sigma_c)
% Gain from data range
p0(5) = log(max(y_range, 0.1));                  % log(gain)
% Surround weight small and subtractive as in example (y - 0.1*y2)
p0(6) = -0.1;                                     % w_s in (-1, 0]

% Robust Huber loss with adaptive delta
delta = 1.345 * mad(y_norm, 1);
if delta <= 0, delta = 0.1; end

% Optimization bounds
lb = [-inf; log(1e-3); 0; log(1e-3); log(1e-3); -1+1e-6];
ub = [inf; log(10); 10; log(10); log(100); 0];

% Objective function
obj_fun = @(p) objective(p, s_norm, y_norm, delta);

% Optimize with multiple random starts for robustness
n_starts = 6;
best_p = p0; best_fval = inf;

for start = 1:n_starts
    if start > 1
        % Random perturbation for additional starts
        p_init = p0 + 0.2 * randn(size(p0));
        p_init = max(min(p_init, ub), lb);
    else
        p_init = p0;
    end
    
    try
        options = optimset('Display', 'off', 'MaxIter', 1200, 'MaxFunEvals', 6000);
        [p_opt, fval] = fminsearch(@(p) obj_fun(constrain_params(p, lb, ub)), p_init, options);
        p_opt = constrain_params(p_opt, lb, ub);
        
        if fval < best_fval
            best_p = p_opt;
            best_fval = fval;
        end
    catch
        % Continue with other starts if one fails
        continue;
    end
end

% Unpack best parameters and denormalize
[mu_c_norm, sigma_c_norm, mu_s_norm, sigma_s_norm, gain_norm, w_s] = unpack_params(best_p);

% Transform back to original scale (no centering in s normalization)
mu_c = mu_c_norm * s_std;
sigma_c = sigma_c_norm * s_std;
mu_s = mu_s_norm * s_std;
sigma_s = sigma_s_norm * s_std;
gain = gain_norm * y_std;

% Calculate global bias to match data mean
s_mid = mean(s);
center_cdf_mid = normcdf_stable(s_mid, mu_c, sigma_c);
surround_cdf_mid = normcdf_stable(s_mid, mu_s, sigma_s);
model_mid = gain * (center_cdf_mid + w_s * surround_cdf_mid);
global_bias = y_mean - model_mid;

% Generate smooth fitted curve (in original s units)
s_plot = linspace(min(s), max(s), 200)';
y_plot = model_original(s_plot, mu_c, sigma_c, mu_s, sigma_s, gain, global_bias, w_s);

% Pack results
p_struct = struct('mu_c', mu_c, 'sigma_c', sigma_c, 'mu_s', mu_s, 'sigma_s', sigma_s, ...
                  'gain', gain, 'global_bias', global_bias, 'w_s', w_s);

stats = struct('fval', best_fval, 'delta', delta, 'n_starts', n_starts);

    function loss = objective(p, s_n, y_n, delta_huber)
        [mu_c_n, sig_c_n, mu_s_n, sig_s_n, gain_n, w_s_val] = unpack_params(p);
        
        % Compute model prediction (without global bias in normalized space)
        cdf_c = normcdf_stable(s_n, mu_c_n, sig_c_n);
        cdf_s = normcdf_stable(s_n, mu_s_n, sig_s_n);
        y_pred = gain_n * (cdf_c + w_s_val * cdf_s);
        
        % Huber loss
        residuals = y_pred - y_n;
        abs_res = abs(residuals);
        huber_vals = (abs_res <= delta_huber) .* (0.5 * residuals.^2) + ...
                     (abs_res > delta_huber) .* (delta_huber * (abs_res - 0.5 * delta_huber));
        loss = mean(huber_vals);
        
        % Add regularization to prevent extreme parameters
        reg = 1e-6 * sum(p.^2);
        
        loss = loss + reg;
    end

    function [mu_c, sigma_c, mu_s, sigma_s, gain, w_s] = unpack_params(p)
        mu_c = p(1);
        sigma_c = exp(p(2));
        mu_s = mu_c + abs(p(3));  % Ensure mu_s >= mu_c
        sigma_s = sigma_c + exp(p(4));  % Ensure sigma_s >= sigma_c
        gain = exp(p(5));
        w_s = max(min(p(6), 0), -1+1e-6);  % Constrain to (-1, 0]
    end

    function p_constrained = constrain_params(p, lower, upper)
        p_constrained = max(min(p, upper), lower);
    end

    function y_pred = model_original(s_vals, mu_c, sig_c, mu_s, sig_s, gain, global_bias, w_s)
        cdf_c = normcdf_stable(s_vals, mu_c, sig_c);
        cdf_s = normcdf_stable(s_vals, mu_s, sig_s);
        y_pred = global_bias + gain * (cdf_c + w_s * cdf_s);
    end

    function cdf_vals = normcdf_stable(x, mu, sigma)
        % Numerically stable normal CDF
        z = (x - mu) ./ max(sigma, 1e-10);
        cdf_vals = 0.5 * (1 + erf(z / sqrt(2)));
    end
end
