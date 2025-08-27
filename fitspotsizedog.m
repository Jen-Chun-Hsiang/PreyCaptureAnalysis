function out = fitspotsizedog(diam, rate, varargin)
% fitspotsizedog  Fit a Difference-of-Gaussians (DoG) model to RGC
% spot-size tuning (mean firing rate vs spot diameter).
%
%   out = fitspotsizedog(diam, rate)
%   out = fitspotsizedog(..., 'BaselineMode','free'|'subtract', 'UseRobust',true/false, 'DoBootstrap',Nboots, 'SizeRatio',4)
%
% Inputs
%   diam : [N x 1] spot diameters (same units; e.g., microns or deg)
%   rate : [N x 1] mean firing rate (Hz) at each diameter
%
% Name-Value options (all optional)
%   'BaselineMode' : 'free' (default) -> fit baseline 'b'
%                    'subtract'       -> subtract min(rate) and fix b=0
%   'UseRobust'    : true/false (default false). If true, uses Huber-like
%                    weights in residuals (simple iteratively reweighted LS).
%   'DoBootstrap'  : integer (default 0). If >0, bootstrap CIs with that many
%                    resamples (with replacement) over (diam,rate) pairs.
%   'SizeRatio'    : scalar (default 4). Forces sigma_s = SizeRatio * sigma_c
%                    during fitting and inference.
%
% Outputs (struct)
%   out.params : [k_c, sigma_c, k_s, sigma_s, b, mu]
%   out.center_weight_int : k_c * (sigma_c^2) * pi
%   out.surround_weight_int: k_s * (sigma_s^2) * pi
%   out.weight_ratio : (k_s*sigma_s^2) / (k_c*sigma_c^2)
%   out.pred     : fitted curve at input diam
%   out.fun      : function handle R = fun(d) for continuous predictions
%   out.gof      : struct with SSE, RMSE, R2
%   out.bootstrap: struct with bootstrapped parameter percentiles if requested
%   out.options  : copy of options used
%
%   The model includes a shift parameter mu (constrained to [-50, 50]) that 
%   allows the center of the Difference-of-Gaussians to be slightly offset from zero.
%
% Dependencies
%   Optimization Toolbox for lsqcurvefit (preferred).
%   Falls back to fminsearch with parameter transforms if lsqcurvefit not found.

% -------- Parse inputs
p = inputParser;
p.addParameter('BaselineMode', 'free', @(s)ischar(s)||isstring(s));
p.addParameter('UseRobust', false, @(x)islogical(x)||ismember(x,[0,1]));
p.addParameter('DoBootstrap', 0, @(x)isnumeric(x)&&isscalar(x)&&x>=0);
p.addParameter('SizeRatio', 4, @(x)isnumeric(x)&&isscalar(x)&&x>1);
p.parse(varargin{:});
opt = p.Results;

diam = diam(:); rate = rate(:);
assert(numel(diam)==numel(rate) && numel(diam)>=4, 'Need matched diam, rate (>=4 points).');

% Optional baseline handling
b_fix = [];
rate_work = rate;
if strcmpi(opt.BaselineMode,'subtract')
    b0 = min(rate_work);
    rate_work = rate_work - b0;
    b_fix = 0;  % fixed baseline
elseif strcmpi(opt.BaselineMode,'free')
    % fit baseline
    b_fix = [];
else
    error('BaselineMode must be ''free'' or ''subtract''.');
end

% ----- Initial guesses
% Peak response and the diameter at (approx) half peak
[~,~] = max(rate_work);
peak = max(rate_work);
if peak <= 0
    peak = max(rate_work) - min(rate_work);
end
% Rough width guess from where response reaches ~half peak
half_level = min(rate_work) + 0.5*(max(rate_work)-min(rate_work));
[~,idx_half] = min(abs(rate_work - half_level));
d_half = max(diam(idx_half), eps);
sig_c0 = max(d_half/2.5,  eps);  % rule-of-thumb

kc0 = max(peak, 1e-2);
ks0 = 0.5*kc0;

if isempty(b_fix)
    b0 = min(rate);
else
    b0 = b_fix;
end

% Parameter vector is now [kc, sigma_c, ks, b, mu] (sigma_s calculated as SizeRatio * sigma_c)
mu0 = 0;  % Initial guess for shift parameter
x0 = [kc0, sig_c0, ks0, b0, mu0];

% Bounds: kc>=0, ks>=0, sigma_c>0, mu constrained to [-50, 50]
lb = [0,  max(sig_c0*1e-3, 1e-6), 0, -Inf, -50];
ub = [Inf, Inf,                 Inf, Inf,  50];

% Model function - sigma_s is calculated as SizeRatio * sigma_c, with shift parameter mu
model = @(x,dd) x(4) + x(1) .* (1 - exp( - ( ((dd-x(5))/2).^2 ) ./ (2*x(2)^2) )) ...
                     - x(3) .* (1 - exp( - ( ((dd-x(5))/2).^2 ) ./ (2*(opt.SizeRatio*x(2))^2) ));

% Objective with optional robust weights (no longer need sigma_s>sigma_c constraint)
    function res = resid_fun(x, dd, yy)
        yhat = model(x, dd);
        res_core = yy - yhat;
        if opt.UseRobust
            % Huber-like weights with scale = 1.345 * MAD
            s = 1.345 * mad(res_core,1);
            s = max(s, 1e-6);
            w = min(1, abs(res_core)/s);
            res_core = w .* res_core;
        end
        res = res_core;  % No penalty needed since sigma_s is automatically > sigma_c
    end

% ---------- Fit
use_lsq = exist('lsqcurvefit','file')==2;

if use_lsq && isempty(b_fix) && ~opt.UseRobust
    % Cleanest path: lsqcurvefit with bounds (no robust, free baseline)
    opts = optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',1e5,'MaxIterations',1e4);
    % Use lsqcurvefit directly on model function
    [xhat, ~, ~, ~] = lsqcurvefit(model, x0, diam, rate, lb, ub, opts);
else
    % Fallback to fminsearch on squared residuals with param transforms to enforce positivity
    % Transform: kc = exp(a), sig_c=exp(b), ks=exp(c), b = free, mu = free
    to_x = @(z) [exp(z(1)), exp(z(2)), exp(z(3)), z(4), z(5)];
    % init in z
    z0 = [log(max(x0(1),1e-6)), log(max(x0(2),1e-6)), log(max(x0(3),1e-6)), x0(4), x0(5)];
    obj = @(z) sum(resid_fun(to_x(z), diam, rate).^2);
    opts = optimset('Display','off','MaxFunEvals',2e5,'MaxIter',2e5);
    [zhat, ~, ~] = fminsearch(obj, z0, opts);
    xhat = to_x(zhat);
    % Constrain mu to [-50, 50]
    xhat(5) = max(-50, min(50, xhat(5)));
end

% If baseline fixed by subtraction, enforce b=0 and refit kc, sigmas, mu only
if ~isempty(b_fix)
    % Refit with b fixed at 0 using lsqcurvefit if available
    model_b0 = @(x,d) 0 + x(1) .* (1 - exp( - ( ((d-x(4))/2).^2 ) ./ (2*x(2)^2) )) ...
                         - x(3) .* (1 - exp( - ( ((d-x(4))/2).^2 ) ./ (2*(opt.SizeRatio*x(2))^2) ));
    x0b = [xhat(1:3), xhat(5)];  % Only kc, sigma_c, ks, mu
    lbb = [0, 1e-6, 0, -50];
    ubb = [Inf, Inf, Inf, 50];
    if use_lsq && ~opt.UseRobust
        optsb = optimoptions('lsqcurvefit','Display','off','MaxFunctionEvaluations',1e5,'MaxIterations',1e4);
        [x1234, ~] = lsqcurvefit(model_b0, x0b, diam, rate_work, lbb, ubb, optsb);
    else
        to_xb = @(z) [exp(z(1)), exp(z(2)), exp(z(3)), z(4)];
        z0b = [log(max(x0b(1),1e-6)), log(max(x0b(2),1e-6)), log(max(x0b(3),1e-6)), x0b(4)];
        objb = @(z) sum( (rate_work - model_b0(to_xb(z), diam)).^2 );
        optsb = optimset('Display','off','MaxFunEvals',2e5,'MaxIter',2e5);
        zhatb = fminsearch(objb, z0b, optsb);
        x1234 = to_xb(zhatb);
        % Constrain mu to [-50, 50]
        x1234(4) = max(-50, min(50, x1234(4)));
    end
    xhat = [x1234(1:3), 0, x1234(4)];  % [kc, sigma_c, ks, b=0, mu]
end

% Predictions, GOF
pred = model(xhat, diam);
res  = rate - pred;
SSE  = sum(res.^2);
SST  = sum( (rate - mean(rate)).^2 );
R2   = 1 - SSE/max(SST, eps);
RMSE = sqrt(mean(res.^2));

kc = xhat(1); sigc = xhat(2); ks = xhat(3); b = xhat(4); mu = xhat(5);
sigs = opt.SizeRatio * sigc;  % Calculate sigma_s from the fixed ratio
center_int   = kc * (sigc^2) * pi;
surround_int = ks * (sigs^2) * pi;
ratio        = (ks * sigs^2) / max(kc * sigc^2, eps);

out = struct();
out.params = [kc, sigc, ks, sigs, b, mu];  % Return 6-element array including mu
out.center_weight_int = center_int;
out.surround_weight_int = surround_int;
out.weight_ratio = ratio;
out.pred = pred;
out.fun  = @(d) model(xhat, d);
out.gof  = struct('SSE',SSE,'RMSE',RMSE,'R2',R2);
out.options = opt;

% ---------- Optional bootstrap CIs
if opt.DoBootstrap > 0
    Nb = opt.DoBootstrap;
    P   = zeros(Nb, 6);  % Updated to 6 parameters including mu
    WR  = zeros(Nb, 1);
    rng('default');
    idxAll = (1:numel(diam))';
    for bsi = 1:Nb
        ii = idxAll(randi(numel(diam), numel(diam), 1));
        try
            bb = fitspotsizedog(diam(ii), rate(ii), 'BaselineMode', opt.BaselineMode, 'UseRobust', opt.UseRobust, 'DoBootstrap', 0, 'SizeRatio', opt.SizeRatio);
            P(bsi,:)  = bb.params;
            WR(bsi,1) = bb.weight_ratio;
        catch
            P(bsi,:)  = NaN(1,6);  % Updated to 6 parameters
            WR(bsi,1) = NaN;
        end
    end
    pct = [2.5 50 97.5];
    out.bootstrap.params_pct = prctile(P, pct, 1);
    out.bootstrap.weight_ratio_pct = prctile(WR, pct, 1);
end
end
