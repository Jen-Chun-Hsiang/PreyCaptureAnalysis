function [s_fit, y_fit, p_struct, fit_stats] = fit_DoG_half_CDF(s, y, modelDim)
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
%   For 2D-line accumulation (along diameter), we use an erf-based shape
%   proportional to erf((D - mu) / (sqrt(2)*sigma)).
%
% INPUTS:
%   s - stimulus sizes (vector; spot diameters)
%   y - neural responses (vector, same length as s)
%   modelDim - optional model type: '1d' (default) or '2d' (line accumulation)
%
% OUTPUTS:
%   s_fit    - fine-grained stimulus sizes for smooth curve
%   y_fit    - fitted response values at s_fit
%   p_struct - structure with fitted parameters:
%              .mu (shared mean), .sigma_c (center std)
%              .sigma_s (surround std)
%              .w_s (surround weight, negative), .gain, .global_bias
%              .modelDim ('1d' or '2d')
%   fit_stats- structure with fitting statistics (.fval, .exitflag, .output)
%
% CONSTRAINTS:
%   - mu: can be slightly negative to zero (shared mean)
%   - sigma_s > sigma_c (surround broader than center)
%   - w_s: negative (surround suppresses response)
%   - gain: positive

% PENALTIES (soft):
%   - Overshoot penalty: model(s_i) > y_i discouraged (quadratic)
%   - Sigma ratio penalty: encourage sigma_s >= 2 * sigma_c (quadratic if violated)
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

% Default model dimension
if nargin < 3 || isempty(modelDim)
    modelDim = '1d'; % options: '1d' (default) or '2d' (line accumulation)
end

% Define half-CDF helpers
% 1D half-CDF (starts at mu): shape ~ 0.5*erf((D-mu)/(sqrt(2)*sigma))
half1d = @(D, mu, sigma) max(0, normcdf(D, mu, sigma) - 0.5);
% 2D line accumulation (integrate 2D Gaussian along a line from mu):
% proportional to erf((D-mu)/(sqrt(2)*sigma)); same shape as 1D up to scale
half2d_line = @(D, mu, sigma) max(0, erf((D - mu) ./ max(sqrt(2)*sigma, eps)) );

% Define the full DoG half-CDF model
% Parameters: [mu, sigma_c, sigma_s, w_s, gain, global_bias]
if strcmpi(modelDim, '2d')
    model_fun = @(params, x) params(5) * ( ...
        half2d_line(x, params(1), params(2)) + params(4) * half2d_line(x, params(1), params(3)) ) ...
        + params(6);
else
    model_fun = @(params, x) params(5) * ( ...
        half1d(x, params(1), params(2)) + params(4) * half1d(x, params(1), params(3)) ) ...
        + params(6);
end

% Objective function (sum of squared residuals with penalties)
overshoot_lambda = 100;   % penalty weight for model exceeding data
sigma_ratio_lambda = 200; % penalty weight for enforcing sigma_s >= 2*sigma_c

objective = @(params) ...
    sum((y - model_fun(params, s)).^2) + ... % SSE
    overshoot_lambda * sum(max(0, model_fun(params, s) - y).^2) + ... % pointwise overshoot
    sigma_ratio_lambda * max(0, 2*params(2) - params(3)).^2; % encourage sigma_s >= 2*sigma_c

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
lb = [0,         max(s_range,eps)*0.01,  max(s_range,eps)*0.05, -1, 0.0,  -inf];
ub = [max(s),    max(s)*2,               max(s)*10,              0,  inf,  inf];

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
if strcmpi(modelDim, '2d')
    p_struct.modelDim = '2d';
else
    p_struct.modelDim = '1d';
end


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
