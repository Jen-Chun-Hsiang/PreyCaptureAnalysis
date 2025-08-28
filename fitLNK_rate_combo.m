function [params, rate_hat, a_traj, fval] = fitLNK_rate_combo(x, x_s, y_rate, dt, varargin)
% fitLNK_rate_combo  Separate LNK kinetics for x and x_s, combined at readout.
%
% Model:
%   a_{t+1}   = a_t   + dt * (alpha_d   * F(x_t   - theta   ) - a_t  ) / tau
%   a_s{t+1}  = a_s{t}+ dt * (alpha_d_s * F(x_s_t - theta_s ) - a_s{t}) / tau_s
%   den   = sigma0   + alpha   * a
%   den_s = sigma0_s + alpha_s * a_s
%   add   = beta * (a + a_s)
%   ỹ_t   = x_t ./ den_t + w_xs * x_s_t ./ den_s_t + add_t + b_out
%   r_t   = outNL( g_out * ỹ_t )
%
% Inputs
%   x, x_s : (T x 1) sequences
%   y_rate : (T x 1) measured firing rate (Hz)
%   dt     : scalar, seconds per bin
%
% Name/Value options
%   'Init'      : struct with fields tau, alpha_d, sigma0, alpha, theta,
%                 tau_s, alpha_d_s, sigma0_s, alpha_s, theta_s,
%                 beta, b_out, g_out, w_xs
%   'MaxIter'   : optimizer iterations (default 400)
%   'Weights'   : (T x 1) per-sample weights for MSE (default ones(T,1))
%   'Robust'    : 'none' (default) or 'huber'
%   'Delta'     : Huber delta (default 1.0)
%   'OutputNL'  : 'softplus' (default) or 'linear'
%   'Ridge'     : L2 penalty on params (default 0)
%
% Outputs
%   params   : struct of fitted parameters
%   rate_hat : (T x 1) fitted rate (Hz)
%   a_traj   : (T x 2) [a, a_s] kinetic state trajectories
%   fval     : final objective value

% ---- options ----
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
init = defaultInitRate_combo(x, x_s, y_rate, dt);
fn = fieldnames(init);
for i = 1:numel(fn)
    if isfield(opts.Init, fn{i}), init.(fn{i}) = opts.Init.(fn{i}); end
end
p0 = packParams_combo(init);

% ---- objective ----
obj = @(p) objectiveAnalog_combo(p, x, x_s, y_rate, dt, opts);

% ---- optimize ----
haveCON  = exist('fmincon','file')==2;

% Param order:
% 1: log(tau), 2: log(alpha_d), 3: log(sigma0), 4: log(alpha),
% 5: beta, 6: b_out, 7: log(g_out), 8: log(theta+1), 9: w_xs,
% 10: log(tau_s), 11: log(alpha_d_s), 12: log(sigma0_s),
% 13: log(alpha_s), 14: log(theta_s+1)

lb = [ ...
    log(1e-3);  % tau
    log(1e-6);  % alpha_d
    log(1e-9);  % sigma0
    log(1e-9);  % alpha
    -10;        % beta
    -100;        % b_out
    log(1e-9);  % g_out
    log(1);     % theta+1
    -2;         % w_xs (suppressive by default)
    log(1e-3);  % tau_s
    log(1e-6);  % alpha_d_s
    log(1e-9);  % sigma0_s
    log(1e-9);  % alpha_s
    log(1)      % theta_s+1
];
ub = [ ...
    log(10);            % tau
    log(10);            % alpha_d
    log(10);            % sigma0
    log(10);            % alpha
    10;                 % beta
    100;                 % b_out
    log(100);           % g_out
    log(max(x)+1+eps);  % theta+1
    0;                  % w_xs
    log(10);            % tau_s
    log(10);            % alpha_d_s
    log(10);            % sigma0_s
    log(10);            % alpha_s
    log(max(x_s)+1+eps) % theta_s+1
];

if haveCON
    disp('Using fmincon');
    o = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',opts.MaxIter);
    [p_hat,fval] = fmincon(obj, p0, [],[],[],[],lb,ub,[], o);
else
    [p_hat,fval] = fminsearch(obj, p0);
end

% ---- unpack & forward ----
params = unpackParams_combo(p_hat);
[rate_hat, a_mat] = forwardLNK_rate_combo(x, x_s, params, dt, opts.OutputNL);
a_traj = a_mat;

end % main

% ----------------- helpers -----------------

function init = defaultInitRate_combo(x, x_s, y, dt)
ypos = max(y,0);
rx  = mean(abs(x(x~=0)));     if isempty(rx),  rx  = 1; end
rxs = mean(abs(x_s(x_s~=0))); if isempty(rxs), rxs = 1; end
ry  = mean(ypos(ypos>0));     if isempty(ry),  ry  = 1; end

% x pathway
init.tau     = max(rand*0.5, 1e-6);
init.alpha_d = max(rand*0.5, 1e-6);
init.theta   = max(0, prctile(x, 20));
init.sigma0  = 0.5*rand + 0.05*std(x);
init.alpha   = 0.05*rand;

% x_s pathway
init.tau_s     = max(rand*0.5, 1e-6);
init.alpha_d_s = max(rand*0.5, 1e-6);
init.theta_s   = max(0, prctile(x_s, 20));
init.sigma0_s  = 0.5*rand + 0.05*std(x_s);
init.alpha_s   = 0.05*rand;

% shared / readout
init.beta    = -0.2*rand;
init.b_out   = 1 + rand;
init.g_out   = max(ry/(rx+rxs+eps), 0.5);
init.w_xs    = -0.1*rand; % default to small suppressive mix
end

function p = packParams_combo(s)
p = [
    log(max(s.tau,      1e-6));
    log(max(s.alpha_d,  1e-6));
    log(max(s.sigma0,   1e-9));
    log(max(s.alpha,    1e-9));
    s.beta;
    s.b_out;
    log(max(s.g_out,    1e-9));
    log(max(s.theta,    0) + 1);   % pair with unpack softplus
    s.w_xs;
    log(max(s.tau_s,     1e-6));
    log(max(s.alpha_d_s, 1e-6));
    log(max(s.sigma0_s,  1e-9));
    log(max(s.alpha_s,   1e-9));
    log(max(s.theta_s,   0) + 1);
];
end

function s = unpackParams_combo(p)
s.tau      = exp(p(1));
s.alpha_d  = exp(p(2));
s.sigma0   = exp(p(3));
s.alpha    = exp(p(4));
s.beta     = p(5);
s.b_out    = p(6);
s.g_out    = exp(p(7));
s.theta    = softplus(p(8));  % ensure >= 0
s.w_xs     = p(9);

s.tau_s     = exp(p(10));
s.alpha_d_s = exp(p(11));
s.sigma0_s  = exp(p(12));
s.alpha_s   = exp(p(13));
s.theta_s   = softplus(p(14)); % ensure >= 0
end

function [rate, a_mat] = forwardLNK_rate_combo(x, x_s, s, dt, outNL)
% Forward pass with separate kinetics for x and x_s

T = numel(x);
a   = zeros(T,1); % for x
a_s = zeros(T,1); % for x_s

for t = 1:T-1
    drive_x  = max(0, x(t)   - s.theta);
    a(t+1)   = a(t)   + dt * (s.alpha_d   * drive_x - a(t))   / s.tau;
    if a(t+1) < 0, a(t+1) = 0; end

    drive_s  = max(0, x_s(t) - s.theta_s);
    a_s(t+1) = a_s(t) + dt * (s.alpha_d_s * drive_s - a_s(t)) / s.tau_s;
    if a_s(t+1) < 0, a_s(t+1) = 0; end
end

den   = s.sigma0   + s.alpha   * a;   den(den<1e-9)   = 1e-9;
den_s = s.sigma0_s + s.alpha_s * a_s; den_s(den_s<1e-9) = 1e-9;

add = s.beta * (a + a_s);

ytilde = x ./ den + s.w_xs * (x_s ./ den_s) + add + s.b_out;

switch lower(outNL)
    case 'softplus'
        rate = softplus(s.g_out * ytilde);
    case 'linear'
        rate = s.g_out * ytilde;
        rate(rate<0) = 0;
    otherwise
        error('OutputNL must be softplus or linear.');
end

a_mat = [a, a_s];
end

function [f, grad] = objectiveAnalog_combo(p, x, x_s, y, dt, opts)
s = unpackParams_combo(p);
[r, ~] = forwardLNK_rate_combo(x, x_s, s, dt, opts.OutputNL);

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

% small L2 on params
f = f + opts.Ridge * (p(:).' * p(:));
if nargout>1, grad = []; end
end

function y = softplus(z)
    y = log1p(exp(-abs(z))) + max(z,0);
end 