close all; clear; clc;
%%
data_sets = {'e100724', 'f100724', 'a101224', 'b101224', 'c101224',   'd101224', 'e101224',...
             'b101424', 'c101424', 'd101424', 'e101424', 'a101624',   'b101624', 'd101624', 'e101624',...
             'b101924', 'c101924', 'd101924', 'e101924', 'b103124',   'e103124', 'a110424',...
             'c110424',                       'f110424', 'g110424',   'a110924', 'b110924', 'c110924',...
             'a111224'};
cell_type = {'OFF',      'OFF',    'OFF',      'ON',       'OFF',     'ON',      'OFF',...
             'OFF',      'OFF',    'ON',       'ON',       'ON',      'ON',      'ON',       'ON',...
             'ON',       'OFF',    'OFF',      'OFF',      'ON',      'OFF',     'ON',...
             'ON',                             'ON',       'OFF',     'ON',      'OFF',      'OFF',...
             'ON'};
location =  {'Temporal', 'Temporal','Nasal',   'Nasal',    'Nasal',   'Nasal',   'Nasal',...
             'Temporal', 'Temporal','Temporal','Temporal', 'Nasal',   'Nasal',   'Nasal',    'Nasal',...
             'Nasal',    'Nasal',   'Nasal',   'Nasal',    'Temporal','Temporal','Temporal', ...
             'Temporal',                      'Temporal',  'Temporal','Temporal','Temporal', 'Temporal',...
             'Temporal'};
clear Data 
num_set = length(data_sets);
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar';
Cdat = nan(num_set, 2, 2);
Cbas = nan(num_set, 2);
implement_case_id = 6;
for i = 1:num_set
    file_name = sprintf('%s_moving_bar_processed.mat', data_sets{i});
    Data{i} = load(fullfile(folder_name, file_name));
    file_name = sprintf('%s_%d_moving_bar_fitted.mat', data_sets{i}, implement_case_id);
    load(sprintf('./Results/MovingBar/%s', file_name), 'PredictionResults', 'BaselineCorr');
    Cdat(i, :, :) = PredictionResults;
    Cbas(i, :) = BaselineCorr;
end
%%
cell_type_numeric = cellfun(@(x) strcmp(x, 'ON'), cell_type);
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), location);
%%
is_blurry = 0;
Colors = lines(4);
figure; 
subplot(1, 2, 1);hold on
cavg = mean(Cbas(cell_type_numeric==1, :), 2);
plot(cavg, 'k');
cavg = mean(Cdat(cell_type_numeric==1, 1, 1), 2);
plot(cavg, 'Color',Colors(1, :));
if is_blurry
    cavg = mean(Cdat(cell_type_numeric==1, 2, 1), 2);
    plot(cavg, 'Color',Colors(2, :));
    legend({'Repeat reliability', 'Standard', 'Blurry'})
else
    legend({'Repeat reliability', 'Standard'})
end
xlabel('# cell');
ylabel('Corr. coeff.');
ylim([0 1])

title('ON');

subplot(1, 2, 2);hold on
cavg = mean(Cbas(cell_type_numeric==0, :), 2);
plot(cavg, 'k');
cavg = mean(Cdat(cell_type_numeric==0, 1, 1), 2);
plot(cavg, 'Color',Colors(1, :));
if is_blurry
    cavg = mean(Cdat(cell_type_numeric==0, 2, 1), 2);
    plot(cavg, 'Color',Colors(2, :));
    legend({'Repeat reliability', 'Standard', 'Blurry'})
else
    legend({'Repeat reliability', 'Standard'})
end
xlabel('# cell');
ylabel('Corr. coeff.');
ylim([0 1])

title('OFF');
%%
keyboard;
%% ON & OFF comparison
Fz = 100;
disp_direction = 0;
disp_contrast = 0.33;
disp_bar_witdth = [50, 100, 200, 400, 800];
disp_speeds = [500, 1000, 2000, 4000, 8000];
max_t = 459;
ct = (0:max_t-1)/Fz;
Colors = parula(4);
cell_type_numeric = cellfun(@(x) strcmp(x, 'ON'), cell_type);
for i = 1:length(disp_contrast)
    Trace = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
    figure;
    for j = 1:length(disp_bar_witdth)
        for q = 1:length(disp_speeds)
            subplot(length(disp_bar_witdth), length(disp_speeds), (j-1)*length(disp_speeds)+q); hold on

            for k = 1:num_set
                dir_id = find(Data{k}.dim1_moving_direction == disp_direction);
                ctr_id = find(Data{k}.dim2_contrast == disp_contrast);
                bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j));
                sp_id = find(Data{k}.dim4_speeds == disp_speeds(q));
                csig = squeeze(mean(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
                Trace(j, q, k, :) = csig(1:max_t);
                switch lower(cell_type{k})
                    case 'on'
                        plot(ct, csig(1:max_t), 'Color', [247, 224, 12]/255);
                    case 'off'
                        plot(ct, csig(1:max_t), 'Color', 0.4*ones(1, 3));
                end
            end
            plot(ct, squeeze(mean(Trace(j, q, cell_type_numeric==1, :), 3)), 'Color', [245 182 66]/255, 'LineWidth', 2)
            plot(ct, squeeze(mean(Trace(j, q, cell_type_numeric==0, :), 3)), 'Color', 0*ones(1, 3), 'LineWidth', 2)
            ylim([0 250]);
            xlim([ct(1) ct(end)]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)')
        end
    end
    sgtitle(sprintf('Direction: %d  Contrast: %0.2G', disp_direction, 1-disp_contrast));
end

%% within type Comparison
Fz = 100;
disp_direction = 0;
disp_contrast = 0;
disp_bar_witdth = [50, 100, 200, 400, 800];
disp_speeds = [500, 1000, 2000, 4000, 8000];
max_t = 459;
Disp_Type = 'ON';

cell_type_id = strcmpi(Disp_Type, 'ON');

ct = (0:max_t-1)/Fz;
Colors = parula(4);
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), location);
for i = 1:length(disp_contrast)
    Trace = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
    figure;
    for j = 1:length(disp_bar_witdth)
        for q = 1:length(disp_speeds)
            subplot(length(disp_bar_witdth), length(disp_speeds), (j-1)*length(disp_speeds)+q); hold on

            for k = 1:num_set
                dir_id = find(Data{k}.dim1_moving_direction == disp_direction);
                ctr_id = find(Data{k}.dim2_contrast == disp_contrast);
                bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j));
                sp_id = find(Data{k}.dim4_speeds == disp_speeds(q));
                csig = squeeze(mean(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
                Trace(j, q, k, :) = csig(1:max_t);
                if strcmpi(cell_type{k}, Disp_Type)
                    switch lower(location{k})
                        case 'temporal'
                            plot(ct, csig(1:max_t), 'Color', [66 182 245]/255);
                        case 'nasal'
                            plot(ct, csig(1:max_t), 'Color', [247 153 205]/255);
                    end
                end
            end
            plot(ct, squeeze(mean(Trace(j, q, cell_type_numeric==cell_type_id & location_type_numeric == 1, :), 3)), 'Color', [27 59 242]/255, 'LineWidth', 2)
            plot(ct, squeeze(mean(Trace(j, q, cell_type_numeric==cell_type_id & location_type_numeric == 0, :), 3)), 'Color', [242 27 145]/255, 'LineWidth', 2)
            ylim([0 250]);
            xlim([ct(1) ct(end)]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)')
        end
    end
    sgtitle(sprintf('%s Direction: %d  Contrast: %0.2G', Disp_Type, disp_direction, 1-disp_contrast));
end
