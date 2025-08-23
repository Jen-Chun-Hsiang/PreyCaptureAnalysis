function [s_fit, y_fit, p_struct, fit_stats] = fit_DoG_half_CDF(s, y)
% FIT_DOG_HALF_CDF Fit a Difference-of-Gaussians model using half-CDF functions
%
% DESCRIPTION:
%   Fits a model where the response is based on the cumulative distribution 
%   function (CDF) of two Gaussians, but only considering the "upper half" 
%   (starting from the mean). This effectively discards the influence of 
%   stimuli smaller than each Gaussian's mean.
%
% MODEL:
%   y = gain * (half_cdf(s, mu, sigma_c) + w_s * half_cdf(s, mu, sigma_s)) + global_bias
%   
%   where half_cdf(s, mu, sigma) = max(0, normcdf(s, mu, sigma) - 0.5)
%
% INPUTS:
%   s - stimulus sizes (vector)
%   y - neural responses (vector, same length as s)
%
% OUTPUTS:
%   s_fit    - fine-grained stimulus sizes for smooth curve
%   y_fit    - fitted response values at s_fit
%   p_struct - structure with fitted parameters:
%              .mu (shared mean), .sigma_c (center std)
%              .sigma_s (surround std)
%              .w_s (surround weight, negative), .gain, .global_bias
%   fit_stats- structure with fitting statistics (.fval, .exitflag, .output)
%
% CONSTRAINTS:
%   - mu: can be slightly negative to zero (shared mean)
%   - sigma_s > sigma_c (surround broader than center)
%   - w_s: negative (surround suppresses response)
%   - gain: positive
%
% AUTHOR: Generated for PreyCaptureRGC analysis
% DATE: August 2025

% Validate inputs
if length(s) ~= length(y)
    error('Input vectors s and y must have the same length');
end

s = s(:);  % ensure column vector
y = y(:);  % ensure column vector

% Remove any NaN or infinite values
valid_idx = isfinite(s) & isfinite(y);
s = s(valid_idx);
y = y(valid_idx);

if length(s) < 5
    error('Need at least 5 valid data points to fit 6 parameters');
end

% Define the half-CDF function
half_cdf = @(x, mu, sigma) max(0, normcdf(x, mu, sigma) - 0.5);

% Define the full DoG half-CDF model
% Parameters: [mu, sigma_c, sigma_s, w_s, gain, global_bias]
model_fun = @(params, x) params(5) * (...
    arrayfun(@(xi) half_cdf(xi, params(1), params(2)), x) + ...
    params(4) * arrayfun(@(xi) half_cdf(xi, params(1), params(3)), x) ...
    ) + params(6);

% Objective function (sum of squared residuals with penalty)
y_max = max(y);
penalty_weight = 1000;  % Large penalty weight for overshooting

objective = @(params) sum((y - model_fun(params, s)).^2) + ...
    penalty_weight * sum(max(0, model_fun(params, s) - y_max).^2);

% Initial parameter estimates
s_range = max(s) - min(s);

% Initial guesses
mu_init = 0;                       % shared mean at 0
sigma_c_init = s_range * 0.1;      % center std (narrow)
sigma_s_init = s_range * 0.8;      % surround std (broader)
w_s_init = -0.3;                   % surround weight (negative)
gain_init = 100;                   % initial gain estimate
bias_init = 0;                     % baseline offset

p0 = [mu_init, sigma_c_init, sigma_s_init, w_s_init, gain_init, bias_init];

% Parameter bounds
% [mu, sigma_c, sigma_s, w_s, gain, global_bias]
lb = [0,             s_range*0,   s_range*0, -1, 0.01, -inf];
ub = [s_range*0.00,  s_range,   s_range*10,   0,  inf,  inf];

% Nonlinear constraints: sigma_s > sigma_c
nonlcon = @(params) deal([], params(2) - params(3));  % sigma_c - sigma_s <= 0

% Optimization options
options = optimoptions('fmincon', ...
    'Display', 'off', ...
    'MaxIterations', 1000, ...
    'MaxFunctionEvaluations', 5000, ...
    'TolFun', 1e-8, ...
    'TolX', 1e-8);

% Perform optimization
try
    [p_opt, fval, exitflag, output] = fmincon(objective, p0, [], [], [], [], ...
        lb, ub, nonlcon, options);
catch ME
    warning('fit_DoG_half_CDF:OptimizationFailed', 'Optimization failed: %s. Using initial guess.', ME.message);
    p_opt = p0;
    fval = objective(p0);
    exitflag = -1;
    output = struct('message', ME.message);
end

% Extract fitted parameters
p_struct = struct();
p_struct.mu = p_opt(1);
p_struct.sigma_c = p_opt(2);
p_struct.sigma_s = p_opt(3);
p_struct.w_s = p_opt(4);
p_struct.gain = p_opt(5);
p_struct.global_bias = p_opt(6);


% Generate smooth fitted curve, x starts at 0
s_min = 0;
s_max = max(s) + s_range * 0.1;
s_fit = linspace(s_min, s_max, 200)';
y_fit = model_fun(p_opt, s_fit);

% Fit statistics
fit_stats = struct();
fit_stats.fval = fval;
fit_stats.exitflag = exitflag;
fit_stats.output = output;

% Calculate R-squared
y_pred = model_fun(p_opt, s);
ss_res = sum((y - y_pred).^2);
ss_tot = sum((y - mean(y)).^2);
fit_stats.r_squared = 1 - ss_res / ss_tot;

% Calculate root mean square error
fit_stats.rmse = sqrt(mean((y - y_pred).^2));

end
