num_data_per_bin = 200;
fit_type = 'cdf'; % 'sigmoid' or 'cdf' to toggle fit function
x_sample_range = -4:0.1:4; % Default range for sampling x
is_wait_to_show = 1;
close all
figure(1); clf; hold on
NL_params = nan(num_set, 6);
NL_curves = nan(num_set, length(x_sample_range));
for j = 1:num_set
    cla; hold on; % Clear axes and hold for each dataset
    file_name = sprintf('%s.mat', data_sets{j});
    load(fullfile(folder_name, file_name), 'PBs', 'FRs');

    [~, sids] = sort(PBs);
    num_bin = ceil(length(PBs)/num_data_per_bin);
    bin_PBs = nan(num_bin, 1);
    bin_FRs = nan(num_bin, 1);
    for i = 1:num_bin
        if i < num_bin
            cids = ((i-1)*num_data_per_bin+1):i*num_data_per_bin;
        else
            cids = ((i-1)*num_data_per_bin+1):length(PBs);
        end
        bin_PBs(i) = mean(PBs(sids(cids)));
        bin_FRs(i) = mean(FRs(sids(cids)));
    end

    % figure; scatter(bin_PBs, bin_FRs, 15, 'k', 'filled');
    %%
        
    x = bin_PBs;
    y = bin_FRs*Fz;

    x = x./std(x);

     % Toggle between fit functions
    switch fit_type
        case 'cdf'
            fit_func = @(p, x) cdf_norm_scaled(x, p(1), p(2), p(3), p(4));
            p0 = [50;    0.1;    1;  1.1;];
            lb = [-Inf, -Inf,    0,    0];   % Lower bounds: sigma >= 0, offset >= 0
            ub = [Inf,   Inf,  Inf,  Inf]; 
        case 'sigmoid'
            fit_func = @(p, x) scaledSigmoid(x, p(1), p(2), p(3), p(4));
            p0 = [30; -1; 1; 0.1];
        otherwise
            error('Unknown fit_type');
    end

    % Perform the fitting
    [pfit, ~] = lsqcurvefit(fit_func, p0, x(:)', y(:)', lb, ub);
    assert(pfit(3) >= 0, 'Fitting failed: sigma < 0');
    assert(pfit(4) >= 0, 'Fitting failed: offset < 0');

    % Sample x from the provided range
    x_sample = x_sample_range;
    y_fit_sample = fit_func(pfit, x_sample);

    % Estimate center of curve
    switch fit_type
        case 'cdf'
            % For cdf_norm_scaled: center at x = B (pfit(4)), y = A*normcdf((0-mu)/1,mu,1)+B
            center_x = pfit(2);
            center_y = cdf_norm_scaled(center_x, pfit(1), pfit(2), pfit(3), pfit(4));
        case 'sigmoid'
            % For scaledSigmoid: center at x = x0 (pfit(3)), y = L/2 + y0
            center_x = pfit(3);
            center_y = pfit(1)/2 + pfit(4);
    end
    NL_params(j, :) = [pfit(:)' center_x center_y];
    NL_curves(j, :) = y_fit_sample;
    % Plot
    scatter(x, y, 5, 'k', 'filled');
    plot(x_sample, y_fit_sample, 'r', 'LineWidth', 2);
    scatter(center_x, center_y, 50, 'b', 'filled'); % Mark center
    xlabel('generator signal');
    ylabel('Firing rate');
    ylim([0 max(y)]);
    title(sprintf('%d, Fit type: %s, Center x: %.3f, Center y: %.3f\n', j, fit_type, center_x, center_y));
    legend('Data', 'Fit', 'Center');
    %%
    fprintf('Fit type: %s, Center x: %.3f, Center y: %.3f\n', fit_type, center_x, center_y);
    if is_wait_to_show
        pause(2)
    end
end
