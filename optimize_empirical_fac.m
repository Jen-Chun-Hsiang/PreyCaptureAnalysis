function [empirical_fac_opt, opt_diff] = optimize_empirical_fac(rgc_time_std, equivalent_noise_level, nRepeats, r_target, T, initial_fac)
    % Uses fzero to solve r_s(fac) = r_target

    opts = optimset('Display','off');
    empirical_fac_opt = fzero( ...
        @(fac) simulate_r_minus_target(fac, rgc_time_std, equivalent_noise_level, nRepeats, T, r_target), ...
        initial_fac, opts);
    opt_diff = simulate_r_minus_target(empirical_fac_opt, rgc_time_std, equivalent_noise_level, nRepeats, T, r_target);
end

function diff = simulate_r_minus_target(fac, rgc_time_std, equivalent_noise_level, nRepeats, T, r_target)
    % Generates repeats, adds noise = equivalent_noise_level*fac,
    % computes mean pairwise Pearson's r_s, and returns r_s - r_target.

    base_trace = rgc_time_std * randn(1, T);
    traces     = repmat(base_trace, nRepeats*10, 1);
    traces     = traces + equivalent_noise_level * fac * randn(size(traces));
    traces     = max(traces, 0);
    R          = corr(traces');
    mask       = triu(true(size(R)), 1);
    r_s        = mean(R(mask));
    diff       = r_s - r_target;
end
