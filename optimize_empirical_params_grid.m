function [fac_opt, offset_opt] = optimize_empirical_params_grid(rgc_time_std, equivalent_noise_level, nRepeats, r_target, mean_rt, rect_thr, T, fac_range, offset_range)
%OPTIMIZE_EMPIRICAL_PARAMS_GRID  Grid search for fac and offset
%   Finds fac (noise scale) and offset that minimize combined error in
%   repeat-reliability and post-rectification mean.
%
%  [fac_opt, offset_opt] = optimize_empirical_params_grid(
%       rgc_time_std,           % STD of smoothed, rectified rgc_time fluctuations
%       equivalent_noise_level, % pre-rectification noise STD x
%       nRepeats,               % number of repeats to simulate per parameter set
%       r_target,               % desired repeat-reliability
%       mean_rt,                % observed mean of smoothed, rectified rgc_time
%       rect_thr,               % rectification threshold used
%       T,                      % number of timepoints per trace
%       fac_range,              % vector of fac values to test
%       offset_range            % vector of offset values to test
%  )

    % Pre-allocate error matrix
    nFac = numel(fac_range);
    nOff = numel(offset_range);
    err_mat = zeros(nFac, nOff);

    % Number of simulated traces per evaluation
    nTr = nRepeats * 10;

    % Loop over grid
    for i = 1:nFac
        fac = fac_range(i);
        for j = 1:nOff
            offset = offset_range(j);
            % 1) Generate raw traces around offset
            raw = offset + rgc_time_std * randn(nTr, T);
            % 2) Add noise before rectification
            noisy = raw + equivalent_noise_level * fac * randn(size(raw));
            % 3) Rectify at threshold
            rectified = max(noisy, rect_thr);
            % 4) Compute mean and reliability
            mean_sim = mean(rectified(:));
            R = corr(rectified');
            mask = triu(true(size(R)), 1);
            r_s = mean(R(mask));
            % 5) Combined squared error
            err = (r_s - r_target)^2 + (mean_sim - mean_rt)^2;
            err_mat(i,j) = err;
        end
    end

    % Find minimum error
    [~, idx] = min(err_mat(:));
    [i_opt, j_opt] = ind2sub(size(err_mat), idx);
    fac_opt = fac_range(i_opt);
    offset_opt = offset_range(j_opt);
end
