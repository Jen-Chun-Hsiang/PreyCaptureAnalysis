function [params, rate_hat, a_traj, fval] = fitLNK_rate_two(x, x_s, y_rate, dt, varargin)
% fitLNK_rate_with_xs  Fit a 1-state LNK model with an additional input x_s.
%
% Model (discrete time, single kinetic state):
%   a_{t+1} = a_t + dt * (alpha_d * F(x_t) - a_t) / tau,  with F(x) = max(0, x - theta)
%   ỹ_t     = x_t + w_xs * x_s_t + [other terms as in LNK]
%   r_t     = outNL( g_out * ỹ_t )   (default outNL = softplus)
%
% INPUTS
%   x      : (T x 1) generator/drive
%   x_s    : (T x 1) secondary input
%   y_rate : (T x 1) measured firing rate (Hz)
%   dt     : scalar, seconds per bin
%
% NAME/VALUE OPTIONS
%   'Init'      : struct with fields tau, alpha_d, sigma0, alpha, beta, b_out, g_out, theta, w_xs
%   'MaxIter'   : optimizer iterations (default 400)
%   'Weights'   : (T x 1) per-sample weights for MSE (default ones(T,1))
%   'Robust'    : 'none' (default) or 'huber'
%   'Delta'     : Huber delta (default 1.0)
%   'OutputNL'  : 'softplus' (default) or 'linear'
%   'Ridge'     : L2 penalty on params (default 0)
%
% OUTPUTS
%   params   : struct of fitted parameters
%   rate_hat : (T x 1) fitted rate (Hz)
%   a_traj   : (T x 1) fitted kinetic state trajectory
%   fval     : final objective value

opts.Init     = struct();
opts.MaxIter  = 400;
opts.Weights  = [];
opts.Robust   = 'none';
opts.Delta    = 1.0;
opts.OutputNL = 'softplus';
opts.Ridge    = 0;
if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end
end

x = x(:); x_s = x_s(:); y_rate = y_rate(:);
T = numel(x);
if isempty(opts.Weights), opts.Weights = ones(T,1); else, opts.Weights = opts.Weights(:); end

% ---- init ----
init = defaultInitRate_xs(x, x_s, y_rate, dt);
fn = fieldnames(init);
for i = 1:numel(fn)
    if isfield(opts.Init, fn{i}), init.(fn{i}) = opts.Init.(fn{i}); end
end
p0 = packParams_xs(init);

% ---- objective ----
obj = @(p) objectiveAnalog_xs(p, x, x_s, y_rate, dt, opts);

haveCON  = exist('fmincon','file')==2;

lb = [log(1e-3); log(1e-6); log(1e-6); log(1e-6); -10; -10; log(1e-6); log(1); -2];
ub = [log(10);   log(10);   log(10);   log(10);   10;   10;   log(100);  log(max(x)+1); 0];

if haveCON
    o = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',opts.MaxIter);
    [p_hat,fval] = fmincon(obj, p0, [],[],[],[],lb,ub,[], o);
else
    [p_hat,fval] = fminsearch(obj, p0);
end

% ---- unpack & forward ----
params = unpackParams_xs(p_hat);
[rate_hat, a_traj] = forwardLNK_rate_xs(x, x_s, params, dt, opts.OutputNL);

end % main

% ----------------- helpers -----------------

function init = defaultInitRate_xs(x, x_s, y, dt)
ypos = max(y,0);
rx = mean(abs(x(x~=0)));  if isempty(rx), rx = 1; end
rxs = mean(abs(x_s(x_s~=0))); if isempty(rxs), rxs = 1; end
ry = mean(ypos(ypos>0));  if isempty(ry), ry = 1; end
init.tau     = rand;
init.alpha_d = rand;
init.theta   = max(0, prctile(x, 20));
init.sigma0  = 0.1 + 0.05*std(x);
init.alpha   = 0.01;
init.beta    = -0.1;
init.b_out   = 1;
init.g_out   = max(ry/(rx+rxs+eps), 0.5);
init.w_xs    = 1;
end

function p = packParams_xs(s)
p = [
    log(max(s.tau,    1e-6));
    log(max(s.alpha_d,1e-6));
    log(max(s.sigma0, 1e-9));
    log(max(s.alpha,  1e-9));
    s.beta;
    s.b_out;
    log(max(s.g_out,  1e-9));
    log(max(s.theta,  0) + 1);
    s.w_xs;
];
end

function s = unpackParams_xs(p)
s.tau     = exp(p(1));
s.alpha_d = exp(p(2));
s.sigma0  = exp(p(3));
s.alpha   = exp(p(4));
s.beta    = p(5);
s.b_out   = p(6);
s.g_out   = exp(p(7));
s.theta   = softplus(p(8));
s.w_xs    = p(9);
end

function [rate, a] = forwardLNK_rate_xs(x, x_s, s, dt, outNL)
T = numel(x);
a = zeros(T,1);
for t = 1:T-1
    drive = max(0, x(t) - s.theta);
    a(t+1) = a(t) + dt * (s.alpha_d * drive - a(t)) / s.tau;
    if a(t+1) < 0, a(t+1) = 0; end
end
den = s.sigma0 + s.alpha * a;  den(den<1e-9) = 1e-9;
add = s.beta * a;
ytilde = x ./ den + s.w_xs * x_s ./ den + add + s.b_out;

switch lower(outNL)
    case 'softplus'
        rate = softplus(s.g_out * ytilde);
    case 'linear'
        rate = s.g_out * ytilde;
        rate(rate<0) = 0;
    otherwise
        error('OutputNL must be softplus or linear.');
end
end

function [f, grad] = objectiveAnalog_xs(p, x, x_s, y, dt, opts)
s = unpackParams_xs(p);
[r, ~] = forwardLNK_rate_xs(x, x_s, s, dt, opts.OutputNL);

res = r - y; w = opts.Weights;
switch lower(opts.Robust)
    case 'none'
        f = mean(w .* (res.^2));
    case 'huber'
        d = opts.Delta;
        a = abs(res);
        hub = (a<=d).* (0.5*res.^2) + (a>d).* (d*(a - 0.5*d));
        f = mean(w .* hub);
    otherwise
        error('Robust must be "none" or "huber".');
end

f = f + opts.Ridge * (p(:).' * p(:));
if nargout>1, grad = []; end
end

function y = softplus(z)
y = log1p(exp(-abs(z))) + max(z,0);
end