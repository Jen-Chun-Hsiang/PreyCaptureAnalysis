num_lnk_test = 10;
lnk_corr_temp = nan(num_lnk_test, 1);
factor_scaling = 1e8;
clear prm_temp
tic
for lnk = 1:num_lnk_test
    fprintf('Fitting LNK model %d/%d... %.2f s\n', lnk, num_lnk_test, toc);
    [prm, r_hat, a_hat, fval] = fitLNK_rate(sim(:)*factor_scaling, exp(:), 0.01, ...
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

clear prm_temp
if is_fitLNK_rate_two
    lnk_corr_temp_s = nan(num_lnk_test, 1);
    tic
    for lnk = 1:num_lnk_test
        fprintf('Fitting LNK model %d/%d... %.2f s\n', lnk, num_lnk_test, toc);
        [prm_s, r_hat_s, a_hat, fval] = fitLNK_rate_two(sim(:)*factor_scaling, sim_s(:)*factor_scaling, exp(:), 0.01, ...
            'OutputNL','softplus', 'Robust','huber', 'Delta',1.0, 'MaxIter',500);
        if lnk == 1
            r_hat_temp_s = nan(num_lnk_test, length(r_hat_s));
        end
        prm_temp{lnk} = prm_s;
        r_hat_temp_s(lnk, :) = r_hat_s;
        lnk_corr_temp_s(lnk) = corr(exp(:), r_hat_s(:));

    end
    [~, bestblk] = max(lnk_corr_temp_s);
    prm_s = prm_temp{bestblk};
    r_hat_s = r_hat_temp_s(bestblk, :);
end

clear prm_temp
if is_fitLNK_divnorm_rate_scw
    lnk_corr_temp_s = nan(num_lnk_test, 1);
    tic
    for lnk = 1:num_lnk_test
        fprintf('Fitting LNK model %d/%d... %.2f s\n', lnk, num_lnk_test, toc);
        [prm_d, r_hat_d, a_hat, fval] = fitLN_divnorm_rate_scw(sim(:)*factor_scaling, sim_s(:)*factor_scaling, exp(:), 0.01, ...
            'OutputNL','softplus', 'Robust','huber', 'Delta',1.0, 'MaxIter',500);
        if lnk == 1
            r_hat_temp_d = nan(num_lnk_test, length(r_hat_d));
        end
        prm_temp{lnk} = prm_d;
        r_hat_temp_d(lnk, :) = r_hat_d;
        lnk_corr_temp_d(lnk) = corr(exp(:), r_hat_d(:));
    end
    [~, bestblk] = max(lnk_corr_temp_d);
    prm_d = prm_temp{bestblk};
    r_hat_d = r_hat_temp_d(bestblk, :);
end

clear prm_temp
if is_fitLNK_rate_scw
    
    lnk_corr_temp_w = nan(num_lnk_test, 1);
    tic
    for lnk = 1:num_lnk_test
        fprintf('Fitting LNK model with CSR constraint %d/%d... %.2f s\n', lnk, num_lnk_test, toc);
        
        % Use CSR_value if it exists, otherwise use default
        if exist('CSR_value', 'var') && ~isempty(CSR_value)
            current_csr = CSR_value;
            if lnk == 1
                fprintf('Using CSR constraint: %.3f\n', current_csr);
            end
        else
            current_csr = 0.09; % default value
            if lnk == 1
                fprintf('Using default CSR constraint: %.3f\n', current_csr);
            end
        end
        
        [prm_w, r_hat_w, a_hat, fval] = fitLNK_rate_scw(sim(:)*factor_scaling, sim_s(:)*factor_scaling, exp(:), 0.01, ...
            'OutputNL','softplus', 'Robust','huber', 'Delta',1.0, 'MaxIter',500, ...
            'CSR', current_csr, 'CSRMetric', 'S_over_C', 'SurroundSign', -1, 'CSRStrength', CSRStrength);
        if lnk == 1
            r_hat_temp_w = nan(num_lnk_test, length(r_hat_w));
        end
        prm_temp{lnk} = prm_w;
        r_hat_temp_w(lnk, :) = r_hat_w;
        lnk_corr_temp_w(lnk) = corr(exp(:), r_hat_w(:));

    end
    [~, bestblk] = max(lnk_corr_temp_w);
    prm_w = prm_temp{bestblk};
    r_hat_w = r_hat_temp_w(bestblk, :);
end
%%
is_save_lnk = 0;
if is_save_lnk
    save_lnk_file_name = sprintf('%s_for_lnk_verification.mat', recording_name);
    save_lnk_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\lnk_verification';
    save(fullfile(save_lnk_folder, save_lnk_file_name), 'sim', 'sim_s', 'exp', ...
        'r_hat', 'prm', 'r_hat_s', 'prm_s', 'r_hat_w', 'prm_w');
end

%%
if is_plot
    x_lim_range = [ct(1) ct(end)];
    close all;
    figure;
    subplot(4, 1, 1)
    plot(ct, exp);
    ylabel('Firing rate (spike/s)');
    yyaxis right
    plot(ct, sim);
    % xlim([0 ct(end)]);
    xlim(x_lim_range);
    xlabel('Time (s)');
    ylabel('linear RF signal (arbi.)');
    title('Data: Experimental firing rate vs Simulation');
    box off;

    subplot(4, 1, 2)
    plot(ct, exp);
    ylabel('Firing rate (spike/s)');
    yyaxis right
    plot(ct, r_hat);
    xlim(x_lim_range);
    xlabel('Time (s)');
    ylabel('Predicted rate (spike/s)');
    title('LNK fit (center only)');
    box off;

    subplot(4, 1, 3)
    plot(ct, exp);
    ylabel('Firing rate (spike/s)');
    yyaxis right
    plot(ct, r_hat_s);
    xlim(x_lim_range);
    xlabel('Time (s)');
    ylabel('Predicted rate (spike/s)');
    title('LNK fit (center + surround)');
    box off;

    subplot(4, 1, 4)
    plot(ct, exp);
    ylabel('Firing rate (spike/s)');
    yyaxis right
    plot(ct, r_hat_w);
    xlim(x_lim_range);
    xlabel('Time (s)');
    ylabel('Predicted rate (spike/s)');
    title('LNK fit (center + surround with CSR constraint)');
    box off;
end