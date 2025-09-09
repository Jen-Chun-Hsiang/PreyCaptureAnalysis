function [params, rate_hat, C_fast, C_slow, fval] = fitLN_divnorm_rate_scw(x, x_s, y_rate, dt, varargin)
% fitLN_divnorm_rate_scw  Fit an LN model with fast+slow divisive contrast normalization.
%
% Model (discrete-time):
%  C_fast_{t+1} = C_fast_t + dt*(|drive_t| - C_fast_t)/tau_cf
%  C_slow_{t+1} = C_slow_t + dt*(|drive_t| - C_slow_t)/tau_cs
%  den_t = sigma0 + gamma_f * C_fast_t + gamma_s * C_slow_t
%  ytilde_t = (drive_t) ./ den_t + b_out
%  r_t = outNL( g_out * ytilde_t )
%  drive_t = x_t + w_xs * x_s_t  (center has implicit weight 1)
%
% INPUTS
%  x, x_s, y_rate, dt : as in fitLNK_rate_scw
%
% NAME/VALUE OPTIONS (defaults in code):
%  'Init' : struct of initial params
%  'MaxIter' : optimizer iterations (default 400)
%  'Weights' : per-sample weights for loss
%  'Robust' : 'none' or 'huber'
%  'Delta' : huber delta
%  'OutputNL' : 'softplus' or 'linear'
%  'Ridge' : L2 penalty on params
%  'MultiStart' : number of random restarts (default 5)
%  'CSR'          : center-surround ratio value (measured independently)
%  'CSRMetric'    : 'S_over_C' (default), 'C_over_S', or 'CSI'
%  'SurroundSign' : -1 (default, inhibitory) or +1 (excitatory)
%  'CSRStrength'  : strength of soft penalty toward CSR-derived w_xs (default 10)
%
% OUTPUTS
%  params : struct of fitted parameters
%  rate_hat : fitted rate (T x 1)
%  C_fast, C_slow : estimated contrast traces (T x 1)
%  fval : final objective value

opts.Init = struct();
opts.MaxIter = 400;
opts.Weights = [];
opts.Robust = 'none';
opts.Delta = 1.0;
opts.OutputNL = 'softplus';
opts.Ridge = 0;
opts.MultiStart = 5; % number of random starts
% --- Center-surround weight options ---
opts.CSR          = [];         % numeric value of your center-surround measure
opts.CSRMetric    = 'S_over_C'; % 'S_over_C', 'C_over_S', or 'CSI'
opts.SurroundSign = -1;         % -1 for inhibitory surround, +1 for excitatory
opts.CSRStrength  = 10;         % strength of soft penalty

% parse varargin
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
init = defaultInit_ln(x, x_s, y_rate, dt);
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

p0 = packParams_ln(init);
obj = @(p) objective_ln(p, x, x_s, y_rate, dt, opts);

haveCON = exist('fmincon','file')==2;

% bounds (in packed p space). Order: log(sigma0), log(gamma_f), log(gamma_s), log(tau_cf), log(tau_cs), b_out, log(g_out), w_xs
lb = [log(1e-9); log(1e-9); log(1e-9); log(1e-3); log(1e-2); -50; log(1e-6); -10];
ub = [log(1e2);  log(1e2);  log(1e2);  log(10);   log(100);   50;  log(100);   10];

if haveCON
    o = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',opts.MaxIter);
else
    o = [];
end

% Multi-start loop
bestf = inf; bestp = p0;
nStarts = max(1, round(opts.MultiStart));
for sidx = 1:nStarts
    if sidx == 1
        pinit = p0;
    else
        pinit = p0 + 0.2 * randn(size(p0));
    end
    try
        if haveCON
            [ptry, ftry] = fmincon(obj, pinit, [],[],[],[], lb, ub, [], o);
        else
            [ptry, ftry] = fminsearch(obj, pinit);
        end
    catch ME
        % fallback: try fminsearch
        try
            [ptry, ftry] = fminsearch(obj, pinit);
        catch
            continue;
        end
    end
    if ftry < bestf
        bestf = ftry; bestp = ptry;
    end
end

p_hat = bestp; fval = bestf;

params = unpackParams_ln(p_hat);
[rate_hat, C_fast, C_slow] = forwardLN_divnorm(x, x_s, params, dt, opts.OutputNL);

% Display final w_xs vs target if CSR was used
if ~isempty(opts.CSR)
    wxs_target = mapCSRtoWxs(opts.CSR, opts.CSRMetric, opts.SurroundSign);
    fprintf('Final w_xs: %.3f (target: %.3f, deviation: %.3f)\n', ...
            params.w_xs, wxs_target, params.w_xs - wxs_target);
end

end

% ---------------- helpers ----------------
function init = defaultInit_ln(x, x_s, y, ~)
% sensible defaults
ypos = max(y,0);
rx = mean(abs(x(x~=0))); if isempty(rx), rx = 1; end
rxs = mean(abs(x_s(x_s~=0))); if isempty(rxs), rxs = 1; end
ry = mean(ypos(ypos>0)); if isempty(ry), ry = 1; end

init.sigma0 = 0.05 * std([x; x_s]) + 1e-3;
init.gamma_f = 0.5 * init.sigma0;    % fast contrast weight
init.gamma_s = 0.1 * init.sigma0;    % slow contrast weight
init.tau_cf = 0.05;  % seconds (fast)
init.tau_cs = 1.0;   % seconds (slow)
init.b_out  = max(0, ry);
init.g_out  = max(ry/(rx + rxs + eps), 0.5);
init.w_xs   = -0.1 * randn();
end

function p = packParams_ln(s)
% pack into vector with logs for positive params
p = [
    log(max(s.sigma0, 1e-9));
    log(max(s.gamma_f, 1e-9));
    log(max(s.gamma_s, 1e-9));
    log(max(s.tau_cf, 1e-6));
    log(max(s.tau_cs, 1e-6));
    s.b_out;
    log(max(s.g_out, 1e-9));
    s.w_xs;
];
end

function s = unpackParams_ln(p)
s.sigma0 = exp(p(1));
s.gamma_f = exp(p(2));
s.gamma_s = exp(p(3));
s.tau_cf = exp(p(4));
s.tau_cs = exp(p(5));
s.b_out = p(6);
s.g_out = exp(p(7));
s.w_xs = p(8);
end

function [rate, C_fast, C_slow] = forwardLN_divnorm(x, x_s, s, dt, outNL)
T = numel(x);
C_fast = zeros(T,1);
C_slow = zeros(T,1);

% init contrast to mean absolute drive
drive = x;
C_fast(1) = mean(abs(drive));
C_slow(1) = mean(abs(drive));

for t = 1:T-1
    C_fast(t+1) = C_fast(t) + dt * (abs(drive(t)) - C_fast(t)) / s.tau_cf;
    C_slow(t+1) = C_slow(t) + dt * (abs(drive(t)) - C_slow(t)) / s.tau_cs;
end

den = s.sigma0 + s.gamma_f * C_fast + s.gamma_s * C_slow;
den(den < 1e-9) = 1e-9;

ytilde = x ./ den + s.w_xs * x_s ./ den + s.b_out;

switch lower(outNL)
    case 'softplus'
        rate = softplus(s.g_out * ytilde);
    case 'linear'
        rate = s.g_out * ytilde;
        rate(rate < 0) = 0;
    otherwise
        error('OutputNL must be softplus or linear.');
end
end

function [f, grad] = objective_ln(p, x, x_s, y, dt, opts)
% compute loss (MSE or Huber) + ridge
s = unpackParams_ln(p);
[r, ~, ~] = forwardLN_divnorm(x, x_s, s, dt, opts.OutputNL);
res = r - y; w = opts.Weights;

switch lower(opts.Robust)
    case 'none'
        f = mean(w .* (res.^2));
    case 'huber'
        d = opts.Delta; a = abs(res);
        hub = (a<=d).* (0.5*res.^2) + (a>d).* (d*(a - 0.5*d));
        f = mean(w .* hub);
    otherwise
        error('Robust must be "none" or "huber".');
end

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
% numerically stable softplus
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
