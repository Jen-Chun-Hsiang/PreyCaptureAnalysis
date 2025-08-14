close all; clear; clc;
%%
set_name = 'latest';  % before081425
is_show_fitted = 0;
switch set_name
     case 'before081425'
        data_sets = {'e100724', 'f100724', 'a101224', 'b101224', 'c101224',   'd101224', 'e101224',...
                    'b101424', 'c101424', 'd101424', 'e101424', 'a101624',   'b101624', 'd101624', 'e101624',...
                    'b101924', 'c101924', 'd101924', 'e101924'};
        cell_type = {'OFF',      'OFF',    'OFF',      'ON',       'OFF',     'ON',      'OFF',...
                    'OFF',      'OFF',    'ON',       'ON',       'ON',      'ON',      'ON',       'ON',...
                    'ON',       'OFF',    'OFF',      'OFF'};
        location =  {'Temporal', 'Temporal','Nasal',   'Nasal',    'Nasal',   'Nasal',   'Nasal',...
                    'Temporal', 'Temporal','Temporal','Temporal', 'Nasal',   'Nasal',   'Nasal',    'Nasal',...
                    'Nasal',    'Nasal',   'Nasal',   'Nasal'};

    case 'latest'
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
end
clear Data 
num_set = length(data_sets);
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar';
Cdat = nan(num_set, 2, 2);
Cbas = nan(num_set, 2);
implement_case_id = 6;
for i = 1:num_set
    file_name = sprintf('%s_moving_bar_processed.mat', data_sets{i});
    Data{i} = load(fullfile(folder_name, file_name));
    if is_show_fitted
        file_name = sprintf('%s_%d_moving_bar_fitted.mat', data_sets{i}, implement_case_id);
        load(sprintf('./Results/MovingBar/%s', file_name), 'PredictionResults', 'BaselineCorr');
        Cdat(i, :, :) = PredictionResults;
        Cbas(i, :) = BaselineCorr;
    end
end
%%
cell_type_numeric = cellfun(@(x) strcmp(x, 'ON'), cell_type);
location_type_numeric = cellfun(@(x) strcmp(x, 'Temporal'), location);
%%
if is_show_fitted
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
end
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
disp_bar_witdth = [50, 100, 200, 400, 800];  % [50, 100, 200, 400, 800]
disp_speeds = [500, 1000, 2000, 4000, 8000]; % [500, 1000, 2000, 4000, 8000]
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
                            plot(ct, csig(1:max_t), 'Color', [180 0 180]/255);
                        case 'nasal'
                            plot(ct, csig(1:max_t), 'Color', [120 0 120]/255);
                    end
                end
            end
            gids_1 = cell_type_numeric==cell_type_id & location_type_numeric == 1;
            gids_2 = cell_type_numeric==cell_type_id & location_type_numeric == 0;
            plot(ct, squeeze(mean(Trace(j, q, gids_1, :), 3)), 'Color', [180 0 180]/255, 'LineWidth', 2)
            plot(ct, squeeze(mean(Trace(j, q, gids_2, :), 3)), 'Color', [120 0 120]/255, 'LineWidth', 2)
            ylim([0 250]);
            xlim([ct(1) ct(end)]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)')
            title(sprintf('%s group 1: n = %d, group 2: n = %d', Disp_Type, sum(gids_1), sum(gids_2)))
        end
    end
    sgtitle(sprintf('%s Direction: %d  Contrast: %0.2G', Disp_Type, disp_direction, 1-disp_contrast));
end

%% within type-location Comparison (4 figures: ON-Temporal, ON-Nasal, OFF-Temporal, OFF-Nasal)
save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Figures\illustrator';
Fz = 100;
disp_direction = 0;
disp_contrast = 0;
disp_bar_witdth = [100];  % [50, 100, 200, 400, 800]
disp_speeds = [500, 8000]; % [500, 1000, 2000, 4000, 8000]
max_t = 459;

ct = (0:max_t-1)/Fz;
group_colors = cat(3, ...
    [180, 0, 180; 120, 0, 120], ... % ON: Temporal, Nasal
    [0, 180, 0; 0, 120, 0]);        % OFF: Temporal, Nasal
group_colors = group_colors / 255;

cell_type_labels = {'ON', 'OFF'};
location_labels = {'Temporal', 'Nasal'};

cell_type_numeric = cellfun(@(x) strcmpi(x, 'ON'), cell_type);
location_type_numeric = cellfun(@(x) strcmpi(x, 'Temporal'), location);

for typeIdx = 1:2
    for locIdx = 1:2
        % Define current group
        Disp_Type = cell_type_labels{typeIdx};
        Disp_Location = location_labels{locIdx};
        % Logical indices for this group
        gids = (cell_type_numeric == (typeIdx==1)) & (location_type_numeric == (locIdx==1));
        
        if sum(gids)==0
            continue; % Skip if no cells in this group
        end

        color_this_group = squeeze(group_colors(locIdx, :, typeIdx));

        for i = 1:length(disp_contrast)
            Trace = nan(length(disp_bar_witdth), length(disp_speeds), num_set, max_t);
            figure('Name', sprintf('%s-%s', Disp_Type, Disp_Location)); % Separate figure for each group
            
            for j = 1:length(disp_bar_witdth)
                for q = 1:length(disp_speeds)
                    subplot(length(disp_bar_witdth), length(disp_speeds), (j-1)*length(disp_speeds)+q); hold on
                    gids_idx = find(gids);
                    for idx = 1:length(gids_idx)
                        k = gids_idx(idx);
                        dir_id = find(Data{k}.dim1_moving_direction == disp_direction);
                        ctr_id = find(Data{k}.dim2_contrast == disp_contrast);
                        bw_id = find(Data{k}.dim3_bar_width == disp_bar_witdth(j));
                        sp_id = find(Data{k}.dim4_speeds == disp_speeds(q));
                        csig = squeeze(mean(Data{k}.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
                        Trace(j, q, k, :) = csig(1:max_t);
                        plot(ct, csig(1:max_t), 'Color', 0.5*ones(1, 3));
                    end
                    % Plot group mean
                    plot(ct, squeeze(mean(Trace(j, q, gids, :), 3)), 'Color', color_this_group, 'LineWidth', 2)
                    ylim([0 250]);
                    yticks(0:50:200);
                    yticklabels({'0', '', '100', '', '200'});
                    switch disp_speeds(q)
                        case 500
                            xlim([ct(1) 3.18]);
                            xticks(0:1:3);  
                            xticklabels({'0',  '1', '2',  '3'});
                        case 8000
                            xlim([ct(1) 1.13]);
                            xticks(0:0.4:1.2);  
                            xticklabels({'0',  '0.4', '0.8',  '1.2'});
                    end
                    
                    xlabel('Time (s)');
                    ylabel('Firing rate (spike/s)')
                    title(sprintf('%s-%s n = %d', Disp_Type, Disp_Location, sum(gids)))
                end
            end
            sgtitle(sprintf('%s-%s Direction: %d  Contrast: %0.2G', Disp_Type, Disp_Location, disp_direction, 1-disp_contrast));
            save_file_name = fullfile(save_folder, sprintf('MovingBar_%s-%s_%d_%0.2G', Disp_Type, Disp_Location, disp_direction, 1-disp_contrast));
            print(gcf, save_file_name, '-depsc', '-vector');
            print(gcf, save_file_name, '-dpng', '-r300');

        end
    end
end