function [empirical_fac_opt, min_err] = optimize_empirical_fac_minimize(rgc_time_std,...
    equivalent_noise_level, nRepeats, r_target, T, max_bin_count, pd2, initial_fac)
%OPTIMIZE_EMPIRICAL_FAC_MINIMIZE  Find fac minimizing reliability error
%
%  [empirical_fac_opt, min_err] = optimize_empirical_fac_minimize(
%       rgc_time_std,           % STD of smoothed, rectified rgc_time fluctuations
%       equivalent_noise_level, % pre-rectification noise STD x
%       nRepeats,               % number of repeats to simulate per evaluation
%       r_target,               % desired repeat-reliability
%       T,                      % number of timepoints per trace
%       initial_fac             % initial guess for fac
%  )
%
% Uses fminsearch to minimize the squared difference (r_s - r_target)^2.

    % Wrapper for fminsearch
    opts = optimset('Display','off', 'TolX',1e-6, 'TolFun',1e-6);
    [empirical_fac_opt, min_err] = fminsearch(@(fac) err_fun(fac, rgc_time_std,...
        equivalent_noise_level, nRepeats, T, max_bin_count, pd2, r_target), initial_fac, opts);
end

function err = err_fun(fac, rgc_time_std, equivalent_noise_level, nRepeats, T, max_bin_count, pd2, r_target)
%ERR_FUN  Compute squared error between simulated r_s and target
    number_gen = 10;
    % Generate base trace and repeats
    base = random(pd2, 1, T);
    base = rgc_time_std*base/std(base);
    % base = rgc_time_std * randn(1, T);
    r_s = nan(number_gen, 1);
    for i =1:number_gen
        traces = repmat(base, nRepeats*5, 1);
        % Add noise before rectification
        traces = traces + equivalent_noise_level * fac * randn(size(traces));
        % Rectify at zero threshold
        traces = max(traces, 0);
        traces = double(traces > mean(traces, 'all'));
        % traces = traces*max_bin_count*40;
        % traces = poissrnd(traces);
        % Compute mean pairwise correlation
        R = corr(traces');
        mask = triu(true(size(R)), 1);
        r_s(i) = mean(R(mask), 'omitnan');
    end
    % Squared error
    m_r_s = mean(r_s)
    err = (m_r_s - r_target)^2;
end