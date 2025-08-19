function [params, rate_hat, a_traj, fval] = fitLNK(x, y, dt, varargin)
% fitLNK  Fit a 1-state LNK (additive + multiplicative) model to an RGC response.
%
% Model (discrete time):
%   a_{t+1} = a_t + dt * (alpha_d * F(x_t) - a_t)/tau,    with F(x)=x (x >= 0)
%   y_tilde = x_t / (sigma0 + alpha * a_t) + beta * a_t + b_out
%   r_t     = softplus(g_out * y_tilde)
%
% INPUTS
%   x  : (T x 1) generator/drive (nonnegative)
%   y  : (T x 1) response; either counts per bin (Poisson) or a firing rate
%   dt : scalar, seconds per bin
%   varargin: name/value pairs:
%       'Loss'    : 'poisson' (default) or 'mse'
%       'Init'    : struct with fields tau, alpha_d, sigma0, alpha, beta, b_out, g_out
%       'MaxIter' : optimizer iterations (default 500)
%
% OUTPUTS
%   params  : struct of fitted parameters (tau, alpha_d, sigma0, alpha, beta, b_out, g_out)
%   rate_hat: (T x 1) fitted rate
%   a_traj  : (T x 1) fitted kinetic state trajectory
%   fval    : final objective value
%
% Dependencies: Optimization Toolbox optional (fminunc/fmincon). Falls back to fminsearch.

opts.Loss    = 'poisson';  % 'poisson' | 'mse'
opts.Init    = struct();
opts.MaxIter = 500;
if ~isempty(varargin)
    for k = 1:2:numel(varargin)
        opts.(varargin{k}) = varargin{k+1};
    end
end

x  = x(:);
y  = y(:);
T  = numel(x);

% ---------- Initialization ----------
init = defaultInit(x, y, dt);
fnames = fieldnames(init);
for i = 1:numel(fnames)
    if isfield(opts.Init, fnames{i})
        init.(fnames{i}) = opts.Init.(fnames{i});
    end
end

% Pack as unconstrained parameters with log transforms for positives
p0 = packParams(init);

% Objective
obj = @(p) objective(p, x, y, dt, opts.Loss);

% ---------- Choose optimizer ----------
useFminunc = exist('fminunc','file') == 2;
useFmincon = exist('fmincon','file') == 2;

if useFminunc
    o = optimoptions('fminunc','Display','off','Algorithm','quasi-newton','MaxIterations',opts.MaxIter);
    [p_hat,fval] = fminunc(obj, p0, o);
elseif useFmincon
    % Unconstrained (we already transform), so empty bounds/constraints
    o = optimoptions('fmincon','Display','off','Algorithm','sqp','MaxIterations',opts.MaxIter);
    [p_hat,fval] = fmincon(obj, p0, [],[],[],[],[],[],[], o);
else
    % Nelder–Mead (no derivatives)
    [p_hat,fval] = fminsearch(obj, p0);
end

% ---------- Unpack & forward pass ----------
params = unpackParams(p_hat);
[rate_hat, a_traj] = forwardLNK(x, params, dt);

end % fitLNK

% ====================== Helpers ======================

function init = defaultInit(x, y, dt)
% crude guesses
ypos = max(y,0);
init.tau     = 0.1;                       % s
init.alpha_d = 1.0;
init.sigma0  = 0.1 + 0.05*std(x);         % baseline divisor
init.alpha   = 0.5;
init.beta    = 0.0;
init.b_out   = 0.0;
% set g_out near ratio of scales
ry   = mean(ypos(ypos>0))+eps; rx = mean(x(x>0))+eps;
init.g_out   = max(ry/rx, 0.5);
end

function p = packParams(s)
% log-transform positive params; others linear
p = [
    log(max(s.tau,    1e-6));
    log(max(s.alpha_d,1e-6));
    log(max(s.sigma0, 1e-9));
    log(max(s.alpha,  1e-9));
    s.beta;
    s.b_out;
    log(max(s.g_out,  1e-9));
];
end

function s = unpackParams(p)
s.tau     = exp(p(1));
s.alpha_d = exp(p(2));
s.sigma0  = exp(p(3));
s.alpha   = exp(p(4));
s.beta    = p(5);
s.b_out   = p(6);
s.g_out   = exp(p(7));
end

function [rate, a] = forwardLNK(x, s, dt)
% One pass through the LNK block, F(x)=x (assumes x>=0)
T  = numel(x);
a  = zeros(T,1);
den = zeros(T,1);
add = zeros(T,1);
for t = 1:T-1
    drive = x(t);                % F(x) = x
    a(t+1) = a(t) + dt * (s.alpha_d * drive - a(t)) / s.tau;
    if a(t+1) < 0, a(t+1) = 0; end
end
den = s.sigma0 + s.alpha * a;
den(den < 1e-9) = 1e-9;
add = s.beta * a;
ytilde = x ./ den + add + s.b_out;
rate   = softplus(s.g_out * ytilde);
end

function [f, grad] = objective(p, x, y, dt, lossType)
% Returns scalar objective; gradient omitted (finite-diff used if needed)
s = unpackParams(p);
[rate, ~] = forwardLNK(x, s, dt);

switch lower(lossType)
    case 'poisson'
        lam = max(rate * dt, 1e-12);
        f = sum(lam - y .* log(lam + 1e-12));      % NLL (up to const)
    case 'mse'
        f = mean((rate - y).^2);
    otherwise
        error('Unknown Loss "%s". Use "poisson" or "mse".', lossType);
end

% Small stabilizer to discourage extreme params (optional)
f = f + 1e-6 * (p(:).' * p(:));
if nargout > 1, grad = []; end
end

function y = softplus(z)
% Numerically-stable softplus
y = log1p(exp(-abs(z))) + max(z,0);
end

% USAGE
%{ 

% Synthetic test
T  = 5000; dt = 0.005;
x  = max(0, randn(T,1)*0.5 + 0.6);                 % nonnegative drive
% Ground truth params
gt.tau=0.08; gt.alpha_d=1.0; gt.sigma0=0.15; gt.alpha=0.8;
gt.beta=0.05; gt.b_out=0.0; gt.g_out=1.2;
[r_true, a_true] = forwardLNK(x, gt, dt);

% Counts from Poisson
y_counts = poissrnd(max(r_true,1e-6)*dt);

% Fit (Poisson)
[prm, r_hat, a_hat, fval] = fitLNK(x, y_counts, dt, 'Loss','poisson', 'MaxIter',600);

% Or fit to analog rate (MSE)
[prm2, r_hat2] = fitLNK(x, r_true, dt, 'Loss','mse', 'MaxIter',600);

fprintf('Fitted tau=%.3f s, alpha_d=%.3f, sigma0=%.3f, alpha=%.3f, beta=%.3f, g_out=%.3f\n',...
    prm.tau, prm.alpha_d, prm.sigma0, prm.alpha, prm.beta, prm.g_out);

% Quick plot
t = (0:T-1)*dt;
subplot(3,1,1); plot(t,x); ylabel('x'); title('Drive');
subplot(3,1,2); plot(t,a_true,'k'); hold on; plot(t,a_hat,'r'); ylabel('a(t)'); legend('true','fit');
subplot(3,1,3); plot(t,r_true,'k'); hold on; plot(t,r_hat,'r'); ylabel('rate'); xlabel('time (s)'); legend('true','fit');


%}
