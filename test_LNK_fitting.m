num_lnk_test = 10;
lnk_corr_temp = nan(num_lnk_test, 1);
clear prm_temp
tic
for lnk = 1:num_lnk_test
    fprintf('Fitting LNK model %d/%d... %.2f s\n', lnk, num_lnk_test, toc);
    [prm, r_hat, a_hat, fval] = fitLNK_rate(sim(:)*1e6, exp(:), 0.01, ...
        'OutputNL','softplus', 'Robust','huber', 'Delta',1.0, 'MaxIter',500);
    if lnk == 1
        r_hat_temp = nan(num_lnk_test, length(r_hat));
    end
    prm_temp{lnk} = prm;
    r_hat_temp(lnk, :) = r_hat;
    lnk_corr_temp(lnk) = corr(exp(:), r_hat(:));

end
[~, bestblk] = max(lnk_corr_temp);
prm = prm_temp{bestblk};
r_hat = r_hat_temp(bestblk, :);

%%
close all
figure;
subplot(2, 1, 1)
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, sim);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
box off;

subplot(2, 1, 2)
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, r_hat);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
box off;