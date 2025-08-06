function [nl_fuc, divider] = getNonlinearFunc(PBs, FRs)

Fz = 100;
num_data_per_bin = 200;
[~, sids] = sort(PBs);
num_bin = round(length(PBs)/num_data_per_bin);
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

figure; scatter(bin_PBs, bin_FRs, 15, 'k', 'filled');

x = bin_PBs;
y = bin_FRs*Fz;
divider = std(x);
x = x./divider;
% Define the fitting function handle
%
fit_func = @(p, x) cdf_norm_scaled(x, p(3), p(2), p(1), p(4));
% Initial guesses for scaling parameters
p0 = [30; 1; 1; 0.1];

% Perform the fitting
[pfit, ~] = lsqcurvefit(fit_func, p0, x(:)', y(:)');

% Evaluate the fitted function
y_fit_cdf = cdf_norm_scaled(x, pfit(3), pfit(2), pfit(1), pfit(4));

nl_fuc = @(x) cdf_norm_scaled(x, pfit(3), pfit(2), pfit(1), pfit(4));
 
fit_func = @(p, x) scaledSigmoid(x, p(1), p(2), p(3), p(4));
p0 = [50; -1; 1; 0.1];
% Perform the fitting
[pfit, ~] = lsqcurvefit(fit_func, p0, x(:)', y(:)');
y_fit_sigmoid = scaledSigmoid(x, pfit(1), pfit(2), pfit(3), pfit(4));

figure; hold on
scatter(x, y, 5, 'k', 'filled');
scatter(x, y_fit_cdf, 5, 'b', 'filled');
scatter(x, y_fit_sigmoid, 5, 'r', 'filled');
xlabel('generator signal');
ylabel('Firing rate');
ylim([0 max(y)]);