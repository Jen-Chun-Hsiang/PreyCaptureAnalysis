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
q_ids = 1; % contrast
bw_ids = 1:5; % barwidth
sp_ids = 1:5; % speed
for k = 1:num_repeat
    for q = q_ids
        iid = randperm(5);
        for i = bw_ids
            for j = sp_ids
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
[repeat_id1, repeat_id2] = randomSplit(num_repeat);
x = mean(Data(dr_id, q_ids, bw_ids, sp_ids, repeat_id1, :), 5);
y = mean(Data(dr_id, q_ids, bw_ids, sp_ids, repeat_id2, :), 5);
rmids = isnan(x) | isnan(y);
x(rmids) = [];
y(rmids) = [];
BaselineCorr = corr(x(:), y(:));
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

test_LNK_fitting
PredictionResults( 3:7) = [corr(exp(:), r_hat(:)) corr(exp(:), r_hat_s(:)),...
                           corr(exp(:), sim(:))   corr(exp(:), sim_s(:)),...
                           corr(exp(:), r_hat_w(:))];
LNK_params = prm;
LNK_params_s = prm_s;
LNK_params_w = prm_w;
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
    PredictionResults(i) = corr(exp(:), y(:));
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
        'BaselineCorr', 'LNK_params', 'LNK_params_s', 'LNK_params_w');
end
all_corr(ii, :) = [PredictionResults(1, :) BaselineCorr(1)];
all_SC(ii) = LNK_params_s.w_xs;
%%
% mean(PredictionResults(:, 1, :), [1 3])
% mean(PredictionResults(:, 2, :), [1 3])
%%
% mean(PredictionResults(:, :, 1), [1 2])
% mean(PredictionResults(:, :, 2), [1 2])
