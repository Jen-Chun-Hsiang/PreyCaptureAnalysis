function res = estimate_rf_sizes(rvals, Rmean, Rsem, bigR, draws)
    % Convenience wrapper
    cen = estimate_center_sigma(rvals, Rmean, Rsem, draws);
    sur = estimate_surround_sigma(rvals, Rmean, Rsem, bigR, draws);

    res = struct();
    res.sigma_c = cen.sigma_c;
    res.ci_c    = cen.ci68;
    res.sigma_s = sur.sigma_s;
    res.ci_s    = sur.ci68;

    % Common diameter conventions
    res.diam63_c = 2*sqrt(2)*res.sigma_c;  % 63%-mass diameter
    res.diam63_s = 2*sqrt(2)*res.sigma_s;
    res.fwhm1d_c = 2.355*res.sigma_c;      % 1-D FWHM equivalent
    res.fwhm1d_s = 2.355*res.sigma_s;

    res.notes = {cen.note, sur.note};
end

function out = estimate_center_sigma(rvals, Rmean, Rsem, draws)
    % Bootstrap half-rise estimate for center sigma:
    %   r50 = radius where ascending curve hits 0.5*peak
    %   sigma_c ≈ r50 / 1.177  (since r50 = sqrt(2 ln 2)*sigma)
    %
    % Uses pchip interpolation for smooth monotonic curves.

    if nargin < 4 || isempty(draws), draws = 2000; end

    rvals = rvals(:)'; Rmean = Rmean(:)'; 
    if ~isempty(Rsem), Rsem = Rsem(:)'; end

    % Peak & ascending limb indices
    [~, i_pk] = max(Rmean);
    idx_up = 1:i_pk;
    r_up = rvals(idx_up);

    % Use pchip for smooth monotonic interpolation (much faster than isotonic)
    R_up = Rmean(idx_up);
    % Ensure monotonically increasing using cummax (simple and fast)
    R_up_smooth = cummax(R_up);

    Rpk = max(R_up_smooth);
    target = 0.5*Rpk;

    % Single-shot r50 from smoothed curve
    r50_base = half_rise_radius(r_up, R_up_smooth, target);

    % Bootstrap (optional)
    sigmas = nan(draws,1);
    for t = 1:draws
        if isempty(Rsem)
            R_samp = R_up_smooth;  % no noise
        else
            R_samp = R_up + Rsem(idx_up).*randn(size(R_up));
            % Simple monotonic enforcement using cummax (much faster than isotonic_inc)
            R_samp = cummax(R_samp);
        end
        Rpk_t = max(R_samp);
        r50_t = half_rise_radius(r_up, R_samp, 0.5*Rpk_t);
        sigmas(t) = r50_t / 1.177;
    end

    med = median(sigmas, 'omitnan');
    lo  = quantile(sigmas, 0.16);
    hi  = quantile(sigmas, 0.84);

    out = struct('sigma_c', med, 'ci68', [lo hi], ...
                'r50_base', r50_base, 'note', 'Using cummax for monotonic smoothing');
end

function r50 = half_rise_radius(r_up, R_up, target)
    % Find first radius where R_up crosses 'target' on ascending limb.
    % Assumes r_up is increasing.
    r_up = r_up(:); R_up = R_up(:);
    % Ensure nondecreasing (safety)
    R_up = max(R_up, [R_up(1); cummax(R_up(1:end-1))]);

    ix = find(R_up >= target, 1, 'first');
    if isempty(ix)
        % never reached: extrapolate to last two points
        ix = numel(r_up);
    end
    if ix==1
        r50 = r_up(1);
    else
        r50 = interp1(R_up([ix-1 ix]), r_up([ix-1 ix]), target, 'linear', 'extrap');
    end
end

function out = estimate_surround_sigma(rvals, Rmean, Rsem, bigR, draws)
    % Bootstrap estimator for surround sigma using:
    %   rho = (R1 - R2) / max(Rpk - R1, eps)
    % and solve_sigma_surround(r1,r2,rho)
    %
    % rvals: vector of spot radii (same units as bigR)
    % Rmean: mean response at rvals
    % Rsem:  SEM (or empty -> deterministic sampling)
    % bigR:  vector of large radii to use, e.g., [400 800 1200]
    % draws: number of bootstrap draws (e.g., 2000)
    %
    % Returns struct with fields:
    %   sigma_s, ci68, n_valid, sigmas (all draws), note

    if nargin < 5 || isempty(draws), draws = 2000; end
    if nargin < 4 || isempty(bigR),  bigR  = [400 800 1200]; end

    rvals = rvals(:)'; Rmean = Rmean(:)'; 
    if ~isempty(Rsem), Rsem = Rsem(:)'; end

    % --- Smooth/synthesize peak via quadratic on local 3 points
    [~, i_pk] = max(Rmean);
    nbrs = unique([max(1,i_pk-1), i_pk, min(numel(rvals),i_pk+1)]);
    p = polyfit(rvals(nbrs), Rmean(nbrs), 2);
    if abs(p(1)) > 0
        rpk_star = -p(2)/(2*p(1));
        rpk_star = min(max(rpk_star, min(rvals(nbrs))), max(rvals(nbrs))); % clamp
    else
        rpk_star = rvals(i_pk);
    end
    Rpk_star = polyval(p, rpk_star);

    pairs = [];
    for i=1:numel(bigR)
        for j=1:numel(bigR)
            if bigR(j) > bigR(i), pairs(end+1,:) = [bigR(i) bigR(j)]; end %#ok<AGROW>
        end
    end

    sigmas = nan(draws*numel(pairs),1);
    c = 0; n_valid = 0;

    for t = 1:draws
        % Optionally sample Rpk (here we keep deterministic Rpk_star; change if you have SEM at pk)
        Rpk = Rpk_star;
        for pidx = 1:size(pairs,1)
            r1 = pairs(pidx,1); r2 = pairs(pidx,2);

            R1 = sample_R_at_radius(r1, rvals, Rmean, Rsem);
            R2 = sample_R_at_radius(r2, rvals, Rmean, Rsem);

            if ~(R1 > R2)
                continue;  % need a drop
            end
            if ~(Rpk > R1)
                continue;  % center saturated by r1
            end
            rho = (R1 - R2) / max(Rpk - R1, 1e-12);

            c = c+1;
            sigmas(c) = solve_sigma_surround(r1, r2, rho);
            n_valid = n_valid + 1;
        end
    end
    sigmas = sigmas(1:c);
    if isempty(sigmas)
        out = struct('sigma_s', NaN, 'ci68', [NaN NaN], 'n_valid', 0, ...
                    'sigmas', sigmas, 'note', 'No informative draws; try larger r2.');
        return
    end

    med = median(sigmas);
    lo  = quantile(sigmas, 0.16);
    hi  = quantile(sigmas, 0.84);

    out = struct('sigma_s', med, 'ci68', [lo hi], 'n_valid', n_valid, ...
                'sigmas', sigmas, 'note', '');
end

function R = sample_R_at_radius(r, rvals, Rmean, Rsem)
    % Linear interp to the requested radius; add Gaussian noise if SEM given.
    Rmu = interp1(rvals, Rmean, r, 'linear', 'extrap');
    if isempty(Rsem)
        R = Rmu;
    else
        Rsd = abs(interp1(rvals, Rsem, r, 'linear', 'extrap'));
        R = Rmu + Rsd.*randn();
    end
end

function sigma = solve_sigma_surround(r1, r2, rho, tol, maxit)
    % Solve for surround sigma in:
    % (exp(-r1^2/2s^2)-exp(-r2^2/2s^2)) / (1-exp(-r1^2/2s^2)) = rho
    % Robust bisection with feasibility clipping.
    %
    % Inputs: r1<r2 (same units), rho in [0, r2^2/r1^2 - 1], tol (opt), maxit (opt)
    % Output: sigma (same units as radii)

    if nargin < 4 || isempty(tol),   tol = 1e-6;   end
    if nargin < 5 || isempty(maxit), maxit = 200;  end

    r1 = double(r1); r2 = double(r2); rho = double(rho);
    if r2 <= r1, error('Require r2 > r1.'); end

    rho_max = (r2^2)/(r1^2) - 1.0;
    rho = max(0.0, min(rho, rho_max));  % clip to feasible set

    f = @(s) ((exp(-r1^2./(2*s.^2)) - exp(-r2^2./(2*s.^2))) ./ ...
            (1 - exp(-r1^2./(2*s.^2)) + 1e-12)) - rho;

    lo = r1/6;        % generous lower bound
    hi = 10*r2;       % generous upper bound

    % If hi not high enough (rare), expand
    while f(hi) < 0 && hi < 1e6*r2
        hi = hi*2;
    end

    % Bisection
    for k = 1:maxit
        mid = 0.5*(lo+hi);
        fm  = f(mid);
        if abs(fm) < tol || (hi-lo) < tol*max(1.0, mid)
            sigma = mid; 
            return
        end
        if fm > 0
            hi = mid;
        else
            lo = mid;
        end
    end
    sigma = 0.5*(lo+hi); % graceful fallback
end
