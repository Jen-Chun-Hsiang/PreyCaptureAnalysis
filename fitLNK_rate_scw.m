function [params, rate_hat, a_traj, fval] = fitLNK_rate_scw(x, x_s, y_rate, dt, varargin)
% fitLNK_rate_scw  Fit a 1-state LNK model with soft center-surround weight constraint.
%
% Model (discrete time, single kinetic state):
%   a_{t+1} = a_t + dt * (alpha_d * F(x_t) - a_t) / tau,  with F(x) = max(0, x - theta)
%   den_t   = sigma0 + alpha * a_t
%   ỹ_t     = x_t / den_t + w_xs * x_s_t / den_t + beta * a_t + b_out
%   r_t     = outNL( g_out * ỹ_t )   (default outNL = softplus)
%
% INPUTS
%   x      : (T x 1) generator/drive (center)
%   x_s    : (T x 1) secondary input (surround)
%   y_rate : (T x 1) measured firing rate (Hz)
%   dt     : scalar, seconds per bin
%
% NAME/VALUE OPTIONS
%   'Init'         : struct with fields tau, alpha_d, sigma0, alpha, beta, b_out, g_out, theta, w_xs
%   'MaxIter'      : optimizer iterations (default 400)
%   'Weights'      : (T x 1) per-sample weights for MSE (default ones(T,1))
%   'Robust'       : 'none' (default) or 'huber'
%   'Delta'        : Huber delta (default 1.0)
%   'OutputNL'     : 'softplus' (default) or 'linear'
%   'Ridge'        : L2 penalty on params (default 0)
%   'CSR'          : center-surround ratio value (measured independently)
%   'CSRMetric'    : 'S_over_C' (default), 'C_over_S', or 'CSI'
%   'SurroundSign' : -1 (default, inhibitory) or +1 (excitatory)
%   'CSRStrength'  : strength of soft penalty toward CSR-derived w_xs (default 10)
%
% OUTPUTS
%   params   : struct of fitted parameters
%   rate_hat : (T x 1) fitted rate (Hz)
%   a_traj   : (T x 1) fitted kinetic state trajectory
%   fval     : final objective value

opts.Init         = struct();
opts.MaxIter      = 400;
opts.Weights      = [];
opts.Robust       = 'none';
opts.Delta        = 1.0;
opts.OutputNL     = 'softplus';
opts.Ridge        = 0;
% --- Center-surround weight options ---
opts.CSR          = [];         % numeric value of your center-surround measure
opts.CSRMetric    = 'S_over_C'; % 'S_over_C', 'C_over_S', or 'CSI'
opts.SurroundSign = -1;         % -1 for inhibitory surround, +1 for excitatory
opts.CSRStrength  = 10;         % strength of soft penalty

if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end
end

x = x(:); x_s = x_s(:); y_rate = y_rate(:);
T = numel(x);
if isempty(opts.Weights), opts.Weights = ones(T,1); else, opts.Weights = opts.Weights(:); end

% Validate CSR inputs
if ~isempty(opts.CSR)
    if ~isnumeric(opts.CSR) || numel(opts.CSR) ~= 1
        error('CSR must be a scalar numeric value');
    end
    if ~ismember(lower(opts.CSRMetric), {'s_over_c', 'c_over_s', 'csi'})
        error('CSRMetric must be ''S_over_C'', ''C_over_S'', or ''CSI''');
    end
    if ~ismember(opts.SurroundSign, [-1, 1])
        error('SurroundSign must be -1 (inhibitory) or +1 (excitatory)');
    end
    fprintf('Using CSR constraint: %s = %.3f, strength = %.1f\n', ...
            opts.CSRMetric, opts.CSR, opts.CSRStrength);
end

% ---- init ----
init = defaultInitRate_xs(x, x_s, y_rate, dt);
fn = fieldnames(init);
for i = 1:numel(fn)
    if isfield(opts.Init, fn{i}), init.(fn{i}) = opts.Init.(fn{i}); end
end

% If CSR is provided, initialize w_xs closer to the target
if ~isempty(opts.CSR)
    wxs_target = mapCSRtoWxs(opts.CSR, opts.CSRMetric, opts.SurroundSign);
    init.w_xs = wxs_target + 0.01*randn; % small noise around target
    fprintf('Target w_xs from CSR: %.3f, initialized at: %.3f\n', wxs_target, init.w_xs);
end

p0 = packParams_xs(init);

% ---- objective ----
obj = @(p) objectiveAnalog_scw(p, x, x_s, y_rate, dt, opts);

haveCON = exist('fmincon','file')==2;

lb = [log(1e-3); log(1e-6); log(1e-6); log(1e-6); -10; -10; log(1e-6); log(1); -3.0];
ub = [log(10);   log(10);   log(10);   log(10);   10;   10;   log(100);  log(max(x)+1); 3.0];

if haveCON
    o = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',opts.MaxIter);
    [p_hat,fval] = fmincon(obj, p0, [],[],[],[],lb,ub,[], o);
else
    [p_hat,fval] = fminsearch(obj, p0);
end

% ---- unpack & forward ----
params = unpackParams_xs(p_hat);
[rate_hat, a_traj] = forwardLNK_rate_xs(x, x_s, params, dt, opts.OutputNL);

% Display final w_xs vs target if CSR was used
if ~isempty(opts.CSR)
    wxs_target = mapCSRtoWxs(opts.CSR, opts.CSRMetric, opts.SurroundSign);
    fprintf('Final w_xs: %.3f (target: %.3f, deviation: %.3f)\n', ...
            params.w_xs, wxs_target, params.w_xs - wxs_target);
end

end % main

% ----------------- helpers -----------------

function init = defaultInitRate_xs(x, x_s, y, dt)
ypos = max(y,0);
rx = mean(abs(x(x~=0)));  if isempty(rx), rx = 1; end
rxs = mean(abs(x_s(x_s~=0))); if isempty(rxs), rxs = 1; end
ry = mean(ypos(ypos>0));  if isempty(ry), ry = 1; end
init.tau     = max(rand*0.5, 1e-6);
init.alpha_d = max(rand*0.5, 1e-6);
init.theta   = max(0, prctile(x, 20));
init.sigma0  = 0.5*rand + 0.05*std(x);
init.alpha   = 0.05*rand;
init.beta    = -0.2*rand;
init.b_out   = 1+rand;
init.g_out   = max(ry/(rx+rxs+eps), 0.5);
init.w_xs    = -0.1*rand;
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

function [f, grad] = objectiveAnalog_scw(p, x, x_s, y, dt, opts)
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

% Ridge regularization
f = f + opts.Ridge * (p(:).' * p(:));

% Soft center-surround weight constraint
if ~isempty(opts.CSR)
    wxs_target = mapCSRtoWxs(opts.CSR, opts.CSRMetric, opts.SurroundSign);
    csr_penalty = opts.CSRStrength * (s.w_xs - wxs_target)^2;
    f = f + csr_penalty;
end

if nargout>1, grad = []; end
end

function y = softplus(z)
y = log1p(exp(-abs(z))) + max(z,0);
end

% --- helper for mapping CSR to w_xs ---
function wxs_target = mapCSRtoWxs(CSR, metric, surroundSign)
% Map center-surround ratio to target w_xs value
%
% INPUTS:
%   CSR         : measured center-surround ratio (should be positive magnitude)
%   metric      : 'S_over_C', 'C_over_S', or 'CSI'
%   surroundSign: -1 (inhibitory) or +1 (excitatory)
%
% OUTPUT:
%   wxs_target  : target value for w_xs parameter

% Ensure CSR is treated as a magnitude (positive value)
CSR = abs(CSR);

switch lower(metric)
    case 's_over_c'
        % CSR = |Surround response| / |Center response|
        k = CSR;
    case 'c_over_s'
        % CSR = |Center response| / |Surround response|
        % Convert to S_over_C format
        k = 1 / max(CSR, eps);
    case 'csi'
        % Center-Surround Index: CSI = (C - S) / (C + S)
        % For inhibitory surround, CSI > 0 means C > |S|
        % For excitatory surround, CSI < 0 means S > C
        % Solve for |S|/|C| from CSI
        r = max(min(CSR, 0.999), -0.999); % clamp to avoid division by zero
        k = abs((1 - r) / (1 + r));
    otherwise
        error('Unknown CSRMetric. Use ''S_over_C'', ''C_over_S'', or ''CSI''.');
end

% Apply sign based on surround type
wxs_target = surroundSign * k;
end
