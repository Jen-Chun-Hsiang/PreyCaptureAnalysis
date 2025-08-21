loadFileName = sprintf('%s_moving_bar_processed.mat', response_name);
load(sprintf('./Results/MovingBar/%s', loadFileName), 'Data');

%% save
if is_blurry
    save_file_name = sprintf('%s_%s_moving_bar_simulated_%d_burry%0.3G_%s.mat', recording_name,...
        response_name, implement_case_id, blurry_length, bar_type);
else
    save_file_name = sprintf('%s_%s_moving_bar_simulated_%d_%s.mat', recording_name,...
        response_name, implement_case_id, bar_type);
end
%%
load(sprintf('./Results/MovingBar/%s', save_file_name), 'dim1_moving_direction', 'dim2_contrast',...
    'dim3_bar_width', 'dim4_speeds','dim5_repeats', 'dim6_time', 'resp', 'cntr');
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
ctr = [];
num_repeat = size(Data, 5);
q_ids = 1; % contrast
bw_ids = 1:5; % barwidth
sp_ids = 1:5; % speed
for k = 1:num_repeat
    for q = q_ids
        iid = randperm(5);
        for i = bw_ids
            for j = sp_ids
                csim = squeeze(resp(q, iid(i), j, :));
                ccntr = squeeze(cntr(q, iid(i), j, :));
                cexp = squeeze(Data(dr_id, q, iid(i), j, k, 1:length(csim)));
                snan_ids = find(isnan(csim));
                enan_ids = find(isnan(cexp));
                if ~isempty(enan_ids)
                    if snan_ids(1) > enan_ids(1)
                        csim(snan_ids(1):end) = [];
                        ccntr(snan_ids(1):end) = [];
                        cexp(snan_ids(1):end) = [];
                    else
                        csim(enan_ids(1):end) = [];
                        ccntr(enan_ids(1):end) = [];
                        cexp(enan_ids(1):end) = [];
                    end
                end
                snan_ids = find(isnan(csim));
                if ~isempty(snan_ids)
                    csim(snan_ids) = csim(snan_ids(1)-1);
                    ccntr(snan_ids) = ccntr(snan_ids(1)-1);
                end

                sim = [sim; csim];
                exp = [exp; cexp];
                ctr = [ctr; ccntr];
            end
        end
    end
end
ct = (0:length(sim)-1)/Fz;
contrast_name = 'std';

keyboard;
%% 
[repeat_id1, repeat_id2] = randomSplit(num_repeat);
x = mean(Data(dr_id, q_ids, bw_ids, sp_ids, repeat_id1, :), 5);
y = mean(Data(dr_id, q_ids, bw_ids, sp_ids, repeat_id2, :), 5);
rmids = isnan(x) | isnan(y);
x(rmids) = [];
y(rmids) = [];
BaselineCorr(jj) = corr(x(:), y(:));
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
plot(ct, ctr);
plot(ct, 1./ctr);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('Contrast status (arbi.)');
ylim([0 max(ctr(:))]);
box off;
%%
assert(mean(isnan(sim))<0.01);
assert(mean(isnan(exp))<0.01);
keyboard
test_LNK_fitting
PredictionResults(jj, 3) = corr(exp(:), r_hat(:));
LNK_params = prm;
%%
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
    PredictionResults(jj, i) = corr(exp(:), y(:));
    PredTraces{jj, i, 1} = exp(:);
    PredTraces{jj, i, 2} = y(:);
    figure; hold on
    plot(ct, exp, 'Color', 0.5*ones(1, 3));
    plot(ct, y, 'b');
    xlim([0 ct(end)]);
    xlabel('Time (s)');
    ylabel('Firing rate (spike/s)');
    clc
    fprintf('progress... %d/%d, %d/%d, %d/%d %.2f s\n', ii, num_recording, jj, 2, i, 2, toc);

    %%
    
    save_file_name = sprintf('%s_%d_moving_bar_fitted.mat', recording_name, implement_case_id);
    save(sprintf('./Results/MovingBar/%s', save_file_name), 'PredictionResults', 'PredTraces',...
        'BaselineCorr', 'LNK_params');
end
all_corr(ii, :) = [PredictionResults(1, :) BaselineCorr(1)];
%%
% mean(PredictionResults(:, 1, :), [1 3])
% mean(PredictionResults(:, 2, :), [1 3])
%%
% mean(PredictionResults(:, :, 1), [1 2])
% mean(PredictionResults(:, :, 2), [1 2])
