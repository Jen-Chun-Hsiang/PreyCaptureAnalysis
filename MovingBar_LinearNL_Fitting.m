loadFileName = sprintf('%s_moving_bar_processed.mat', response_name);
load(sprintf('./Results/MovingBar/%s', loadFileName), 'Data');

%% save
movingbar_sim_folder = './Results/MovingBar/';
str1 = sprintf('%s_%s_moving_bar_LN_simulated', recording_name,...
    response_name);


files = dir(fullfile(movingbar_sim_folder, '*.mat'));
matchIdx = find(contains({files.name}, str1) & contains({files.name}, bar_type));

if ~isempty(files)
    if length(matchIdx) > 1
        error('Multiple matching files found.');
    else
        matchedFile = files(matchIdx).name;
        fprintf('Using matched file: %s\n', matchedFile);
    end
else
    error('No matching files found.');
end

%%
load(fullfile(movingbar_sim_folder, matchedFile), 'dim2_contrast',...
    'dim3_bar_width', 'dim4_speeds','dim5_repeats', 'dim6_time', 'resp', 'resp_s');

%%
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_' stim_wn_id '_Retina_1_MovingNoise_1.mat'], 'OLED');
pix2um = OLED.pixelSize;
blurry_length = 50/pix2um; % 0.5*
%% stimulus calculate of contrast
dr_id = 1;
% q = 1; %contra_id
% i = 3; % bar_width_id
% j = 1; % speed_id
ctrs = [1 2/3 1/3]*2;
sim = [];
exp = [];
sim_s = [];
num_repeat = size(Data, 5);
q_ids = 1:1; % contrast
bw_ids = 1:5; % barwidth
sp_ids = 1:5; % speed

trial_table = nan(num_repeat*length(q_ids)*length(bw_ids)*length(sp_ids), 5);
tid = 1;
rpids = [];
trail = [];
for k = 1:num_repeat
    for q = q_ids
        iid = randperm(length(bw_ids));
        for i = bw_ids
            for j = sp_ids
                trial_table(tid, :) = [q, iid(i), j, k, i];

                
                csim = squeeze(resp(q, iid(i), j, :));
                csim_s = squeeze(resp_s(q, iid(i), j, :));
                cexp = squeeze(Data(dr_id, q, iid(i), j, k, 1:length(csim)));
                snan_ids = find(isnan(csim));
                enan_ids = find(isnan(cexp));
                if ~isempty(enan_ids)
                    if snan_ids(1) > enan_ids(1)
                        csim(snan_ids(1):end) = [];
                        csim_s(snan_ids(1):end) = [];
                        cexp(snan_ids(1):end) = [];
                    else
                        csim(enan_ids(1):end) = [];
                        csim_s(enan_ids(1):end) = [];
                        cexp(enan_ids(1):end) = [];
                    end
                end
                snan_ids = find(isnan(csim));
                if ~isempty(snan_ids)
                    csim(snan_ids) = csim(snan_ids(1)-1);
                    csim_s(snan_ids) = csim_s(snan_ids(1)-1);
                end
                trail = [trail; tid*ones(length(csim), 1)];
                rpids = [rpids; k*ones(length(cexp), 1)];
                tid = tid + 1;
                sim = [sim; csim];
                exp = [exp; cexp];
                sim_s = [sim_s; csim_s];
            end
        end
    end
end
ct = (0:length(sim)-1)/Fz;
contrast_name = 'std';

%% 
[repeat_data] = Data(dr_id, q_ids, bw_ids, sp_ids, :, :);
if num_repeat < 2
    BaselineCorr = NaN;
else
    repeat_data = reshape(repeat_data, [], num_repeat, size(repeat_data, 6));
    repeat_data = permute(repeat_data, [2 1 3]);
    repeat_data = reshape(repeat_data, num_repeat, []);
    num_pairs = nchoosek(num_repeat, 2);
    pair_corrs = nan(num_pairs, 1);
    pair_idx = 1;
    for r1 = 1:(num_repeat-1)
        x = repeat_data(r1, :)';
        for r2 = (r1+1):num_repeat
            y = repeat_data(r2, :)';
            valid = ~isnan(x) & ~isnan(y);
            if any(valid)
                pair_corrs(pair_idx) = corr(x(valid), y(valid));
            end
            pair_idx = pair_idx + 1;
        end
    end
    BaselineCorr = mean(pair_corrs, 'omitnan');
end

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

subplot(2, 1, 2);
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
hold on
plot(ct, sim-0.05*sim_s);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('Contrast status (arbi.)');
% ylim([0 max(ctr(:))]);
box off;
%%
assert(mean(isnan(sim))<0.01);
assert(mean(isnan(sim_s))<0.01);
assert(mean(isnan(exp))<0.01);
sim(isnan(sim)) = 0;
sim_s(isnan(sim_s)) = 0;
exp(isnan(exp)) = 0;

is_fitLNK_rate_two = 0;
is_fitLNK_divnorm_rate_scw = 1;
is_fitLNK_rate_scw = 0;
is_fitGainControl = 0;
test_LNK_fitting
[sim_nl, LN_params] = fit_linear_transform(sim*1e6, exp);
% Check if CSR values are available and pass them to the fitting function
if exist('CSR_value', 'var') && exist('CSRStrength', 'var')
    disp('Fitting with CSR constraint...');
    [sim_nl_s, LN_params_s] = fit_linear_transform_with_surround(sim*1e6, exp, sim_s*1e6, 'CSR_value', CSR_value, 'CSRStrength', CSRStrength);
else
    [sim_nl_s, LN_params_s] = fit_linear_transform_with_surround(sim*1e6, exp, sim_s*1e6);
end
assert(length(sim_nl_s) == numel(trail));

if is_fitLNK_rate_two
    r_hat_s_corr = repeatcorr(exp(:), r_hat_s(:), rpids(:));
    LNK_params_s = prm_s;
else
    r_hat_s_corr = nan; 
    LNK_params_s = [];
    LNK_params_s.w_xs = nan;
end
if is_fitLNK_divnorm_rate_scw
    r_hat_d_corr = repeatcorr(exp(:), r_hat_d(:), rpids(:));
    LNK_params_d = prm_d;
else
    r_hat_d_corr = nan; 
    LNK_params_d = [];
    LNK_params_d.w_xs = nan;
end
if is_fitLNK_rate_scw
    r_hat_w_corr = repeatcorr(exp(:), r_hat_w(:), rpids(:));
    LNK_params_w = prm_w;
else
    r_hat_w_corr = nan;
    LNK_params_w = [];
    LNK_params_w.w_xs = nan;
end
PredictionResults( 3:end) = [repeatcorr(exp(:), r_hat(:), rpids(:))   r_hat_s_corr,...
                           repeatcorr(exp(:), sim_nl(:), rpids(:))  repeatcorr(exp(:), sim_s(:), rpids(:)),...
                           r_hat_w_corr  repeatcorr(exp(:), sim(:), rpids(:)),...
                           repeatcorr(exp(:), sim_nl_s(:), rpids(:)), r_hat_d_corr];
LNK_params = prm;



%% Fit gain control model
if is_fitGainControl
    csim = sim;
    csim = csim*1e6;
    gain_params = [200, 100, 0.3, 1, 0.01, 0.01];

    for i = 1:2
        switch i
            case 1
                is_sigmoid = 0;
            case 2
                is_sigmoid = 1;
        end
        tic
        if is_sigmoid
            CostF = @(w) mean((max([gain_control_system_opt_sigmoid([gain_params(1:2) w(1:4)], csim).*w(5)+...
                w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
            [OptW,fval] = fmincon(CostF, ...
                [0.1   1    0.01   0.01  1     0    ], [], [], [], [],...
                [1e-8  0    1e-8   -2    1e-3  -200 ],...
                [1    1e3   2      2     1e7   200  ]);
        else
            CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*w(5)+...
                w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
            [OptW,fval] = fmincon(CostF, ...
                [0.1   1    4   0.01    1     0    ], [], [], [], [],...
                [1e-8  0    1   -4      1e-3  -200 ],...
                [1    1e3   10   4      1e7   200  ]);
        end
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*OptW(5)+OptW(6);
            zeros(1, length(csim))], [], 1);
        ct = (0:length(sim)-1)/Fz;
        PredictionResults(i) = repeatcorr(exp(:), y(:), rpids(:));
        PredTraces{i, 1} = exp(:);
        PredTraces{i, 2} = y(:);
        figure; hold on
        plot(ct, exp, 'Color', 0.5*ones(1, 3));
        plot(ct, y, 'b');
        xlim([0 ct(end)]);
        xlabel('Time (s)');
        ylabel('Firing rate (spike/s)');
        clc
        fprintf('progress... %d/%d, %d/%d %.2f s\n', ii, num_recording, i, 2, toc);

        %%
        save_file_name = sprintf('%s_moving_bar_fitted.mat', recording_name);
        save(sprintf('./Results/MovingBar/%s', save_file_name), 'PredictionResults', 'PredTraces',...
            'BaselineCorr', 'LNK_params', 'LNK_params_s', 'LNK_params_w', 'LN_params', 'LN_params_s', 'LNK_params_d', 'sim_nl_s',...
            'trail', 'trial_table');
    end
else
    % If gain control fitting is disabled, set default values for gain control indices
    PredictionResults(1:2) = NaN;
    PredTraces = cell(2, 2);
    for i = 1:2
        PredTraces{i, 1} = exp(:);
        PredTraces{i, 2} = NaN(size(exp));
    end
    
    save_file_name = sprintf('%s_moving_bar_fitted.mat', recording_name);
    save(sprintf('./Results/MovingBar/%s', save_file_name), 'PredictionResults', 'PredTraces',...
        'BaselineCorr', 'LNK_params', 'LNK_params_s', 'LNK_params_w', 'LN_params', 'LN_params_s', 'LNK_params_d', 'sim_nl_s',...
        'trail', 'trial_table');
end
all_corr(ii, :) = [PredictionResults(1, :) BaselineCorr(1)];
all_SC(ii, :) = [LNK_params_s.w_xs LNK_params_w.w_xs LN_params_s.gamma LNK_params_d.w_xs];
%%
% mean(PredictionResults(:, 1, :), [1 3])
% mean(PredictionResults(:, 2, :), [1 3])
%%
% mean(PredictionResults(:, :, 1), [1 2])
% mean(PredictionResults(:, :, 2), [1 2])


function c = repeatcorr(x, y, repeat_ids)
%REPEATCORR Correlation averaged across repeats identified by repeat_ids.
%   c = REPEATCORR(x, y, repeat_ids) returns the mean Pearson correlation
%   between vectors x and y computed separately for each unique repeat
%   index in repeat_ids. NaN entries are ignored on a per-repeat basis and
%   repeats with fewer than two valid samples are skipped. The returned
%   value is the mean of the per-repeat correlations, omitting NaNs.

    if ~isequal(numel(x), numel(y), numel(repeat_ids))
        error('repeatcorr:InputSizeMismatch', ...
            'Inputs x, y, and repeat_ids must have the same number of elements.');
    end

    x = x(:);
    y = y(:);
    repeat_ids = repeat_ids(:);

    uids = unique(repeat_ids(~isnan(repeat_ids)));
    per_repeat = nan(numel(uids), 1);

    for idx = 1:numel(uids)
        mask = repeat_ids == uids(idx);
        xv = x(mask);
        yv = y(mask);
        valid = ~isnan(xv) & ~isnan(yv);
        if nnz(valid) > 1
            per_repeat(idx) = corr(xv(valid), yv(valid));
        end
    end

    c = mean(per_repeat, 'omitnan');
end
