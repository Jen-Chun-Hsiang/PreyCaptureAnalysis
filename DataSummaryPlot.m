clear close all; clc
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
%% There is a middle point loading (line 234)
keyboard;

%%
% Specify the Excel file path
excelFilePath = 'PreyCaptureRGC_Data.xlsx'; % Update this with your actual file path
resultFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results';
sheetName = 'Recording_Info';

% Read the Excel sheet into a table
dataTable = readtable(excelFilePath, 'FileType', 'spreadsheet', 'Sheet', sheetName);

dataTable.TN_axis_str = dataTable.TN_axis;
TN_values = (strcmp(dataTable.TN_axis, 'T') - 1*strcmp(dataTable.TN_axis, 'N') + 1)* 0.5;
dataTable.TN_axis = TN_values;

dataTable.DV_axis_str = dataTable.DV_axis;
DV_values = (strcmp(dataTable.DV_axis, 'D') - 1*strcmp(dataTable.DV_axis, 'V') + 1)* 0.5;
dataTable.DV_axis = DV_values;

dataTable.Side_str = dataTable.Side;
values = -strcmp(dataTable.Side, 'L') + strcmp(dataTable.Side, 'R');
dataTable.Side = values;

dataTable.Condition_str = dataTable.Condition;
values = strcmp(dataTable.Condition, 'MFA');
dataTable.Condition = values;


subDataTable = dataTable(:, {'Day', 'Side'});
[uniqueRows, ~, ic] = unique(subDataTable, 'rows', 'stable');
dataTable.Retina_id = ic;

%*** check Side reconfiguration

% Verify if "file_name" column exists
if ~ismember('File_name', dataTable.Properties.VariableNames)
    error('The column "file_name" does not exist in the specified sheet.');
end

%%
% Initialize a struct to hold loaded .mat files and other info
loadedData = struct;
DataInfo = nan(height(dataTable), 5);
clear Data_Spot Data_MovStop
Fz = 100;
WinT = [-0.5 0];

% Loop through each row in the table
for rowIndex = 1:height(dataTable)
    clc
    fprintf('loading... %d/%d \n', rowIndex, height(dataTable));
    % Extract the file name from the "file_name" column
    matFileName = dataTable.File_name{rowIndex};

    % for white noise
    loadFileName = fullfile([resultFolder '\MovingWhite'], [matFileName, '.mat']);
    % Verify and load the .mat file
    if exist(loadFileName, 'file')
        DataInfo(rowIndex, :) = [dataTable.Cell_type(rowIndex) dataTable.TN_axis(rowIndex) dataTable.Retina_id(rowIndex),...
            dataTable.Side(rowIndex) dataTable.Condition(rowIndex)];
        loadedFile = load(loadFileName);
        % method 1: Find peaks
        smtstdSTA = imgaussfilt(loadedFile.stdSTA,5);
        smtstdSTA = medfilt2(smtstdSTA, [5, 5]);
        a = smtstdSTA';
        a = 256*(a-min(a(:)))/range(a(:));
        thr = quantile(a(:), 0.99);
        Cent = FastPeakFind(a, thr);
        Cent = reshape(Cent, 2, [])';
        Data_Cent{rowIndex} = Cent;

        % method 2: Gaussian correlation
        smtstdSTA = medfilt2(loadedFile.stdSTA, [3, 3]);
        image = smtstdSTA';
        initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];

        if rowIndex == 1
            data_optimal_params_1 = nan(height(dataTable), length(initial_params));
            data_optimal_params_2 = nan(height(dataTable), length(initial_params));
        end
        mask_1 = ~isnan(loadedFile.masked_stdSTA);
        objective_function = @(params) gaussian_difference(params, image, mask_1);
        optimal_params_1 = fminsearch(objective_function, initial_params);
        data_optimal_params_1(rowIndex, :) = optimal_params_1;

        mask_2 = loadedFile.masked_stdSTA>=0;
        objective_function = @(params) gaussian_difference(params, image, mask_2);
        optimal_params_2 = fminsearch(objective_function, initial_params);
        data_optimal_params_2(rowIndex, :) = optimal_params_2;


        [X, Y] = meshgrid(1:size(image,2), 1:size(image,1));

        if ~exist('Data_GaussCorr', 'var')
            Data_GaussCorr = nan(height(dataTable), 2);
            Data_GaussFit = nan(height(dataTable), numel(optimal_params_1), 2);
        end
        Data_GaussFit(rowIndex, :, 1) = optimal_params_1;
        Data_GaussFit(rowIndex, :, 2) = optimal_params_2;
        gaussian_model = gaussian2d(X, Y, optimal_params_1);
        Data_GaussCorr(rowIndex, 1) = corr(reshape(image(mask_1), [], 1), reshape(gaussian_model(mask_1), [], 1))^2;

        gaussian_model = gaussian2d(X, Y, optimal_params_2);
        Data_GaussCorr(rowIndex, 2) = corr(reshape(image(mask_2), [], 1), reshape(gaussian_model(mask_2), [], 1))^2;

        % get TF
        smtSTAmat = medfilt3(loadedFile.STAmat);
        smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
        tRF = reshape(smtSTAmat, [], size(loadedFile.STAmat, 3))'*loadedFile.masked_stdSTA(:);
        t = WinT(1):1/Fz:WinT(end);
        % figure; plot(t(2:end), tRF, 'k');
        if ~exist('Data_TF', 'var')
            Data_TF = nan(height(dataTable), length(tRF));
        end
        Data_TF(rowIndex, :) = tRF;
    else
        warning('File "%s" does not exist.', matFileName);
    end

    % for moving spot
    loadFileName = fullfile([resultFolder '\MovingSpot'], [matFileName, '_moving_spot_processed.mat']);
    % Verify and load the .mat file
    if exist(loadFileName, 'file')
        loadedFile = load(loadFileName);
        Data_Spot{rowIndex} = loadedFile.Data;
        Data_MovStop{rowIndex} = loadedFile.Moving_end_time_ids;
    else
        warning('File "%s" does not exist.', matFileName);
    end
end
Data_Spot{rowIndex} = [];
matFileName = dataTable.File_name{1};
loadFileName = fullfile([resultFolder '\MovingSpot'], [matFileName, '_moving_spot_processed.mat']);
loadedFile = load(loadFileName);
%%
keyboard;
%%
cids = DataInfo(:, 1)*2+DataInfo(:, 2);
figure; 
sum_params = mean(data_optimal_params_1(cids == 0, :), 1);
Z = gaussian2d(X, Y, sum_params);
imagesc(Z); colorbar;

%%
cent_count = cellfun(@(x) size(x, 1), Data_Cent);
figure;
scatter(cent_count, Data_GaussCorr(:, 1).^2, 15, 'k', 'filled')

%% Gaussian-like RF measurement
figure; hold on
%swarmchart(DataInfo(:, 1), cent_count)
cids = DataInfo(:, 1)*2+DataInfo(:, 2);
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
sids = DataInfo(:, 5) == 0 & (1:size(DataInfo, 1))' ~=32;
swarmchart(cids(sids), Data_GaussCorr(sids, 1).^2, 25, 'k', 'filled');
for i = 1:length(ucids)
    plot(i -1 + [-0.3 0.3], mean(Data_GaussCorr(cids == ucids(i) & sids, 2).^2)*ones(1, 2), 'k');
end
xticks(0:3)
xticklabels({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});
xlabel('Cell category');
ylabel('Gaussian-like index (R2)');
%%
[p,tbl,stats] = kruskalwallis(Data_GaussCorr(sids, 1).^2,cids(sids),'off');
c = multcompare(stats);
%%

figure; 
aids = find(ismember(DataInfo(:, 3), [8 9]));
sids = aids(DataInfo(aids, 5) == 0);
sids([1 end]) = [];
cids = aids(DataInfo(aids, 5) == 1);
plot([Data_GaussCorr(sids, 1).^2 Data_GaussCorr(cids, 1).^2]');
xticks(1:2)
xticklabels({'Control', 'MFA'});
xlim([0.8 2.2])
xlabel('Conditions');
ylabel('Gaussian-like index (R2)');
box off
%% Peak find RF measurement
figure;
%swarmchart(DataInfo(:, 1), cent_count)
swarmchart(DataInfo(:, 1)*2+DataInfo(:, 2), cent_count, 'k', 'filled')
xticks(0:3)
xticklabels({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});
xlabel('Cell category');
ylabel('Number of RF peaks');

%% Peak find RF measurement
clear cent_cells
for i = 1:length(Data_Cent)
    if ~isempty(Data_Cent{i})
        cent_cells(i, :) = mean(Data_Cent{i}, 1);
    end
end

cids = DataInfo(:, 1)*2+DataInfo(:, 2);
ucids = unique(cids);
ucids = ucids(~isnan(ucids));

figure; hold on
for i = ucids(:)'
    switch i
        case 2 % ON-Nasal
            line_color = [92 204 206]/255;
        case 3 % ON-Temporal
            line_color = [255 139 139]/255;
        case 0 % OFF-Nasal
            line_color = [60 91 89]/255;
        case 1 % OFF-Temporal
            line_color = [136 74 75]/255;
    end
    scatter(cent_cells(cids==i & all(~isnan(cent_cells), 2), 1), cent_cells(cids==i & all(~isnan(cent_cells), 2), 2), 25,  line_color, 'filled');
end
plot([0 800], 300*ones(1, 2), '--k');
plot(400*ones(1, 2), [0 600], '--k');
% xticks(0:3);
legend({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});
xlabel('Width (pixel)');
ylabel('Height (pixel)');
xlim([300 500]);
ylim([150 400]);
title('Cell types')
%%
cids = DataInfo(:, 3);
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
colors = lines(length(ucids));
clear retina_cent;
figure; hold on
for i = ucids(:)'
    sids = cids==i & all(~isnan(cent_cells), 2) & (DataInfo(:, 1) == 0) & (DataInfo(:, 5)== 0);
    retina_cent{i} = cent_cells(sids, :);
    scatter(cent_cells(sids, 1), cent_cells(sids, 2),...
        25,  colors(i, :), 'filled');
    % scatter(cent_cells(cids==i & all(~isnan(cent_cells), 2) & (DataInfo(:, 1) == 0 | DataInfo(:, 2) == 0), 1),...
    %     cent_cells(cids==i & all(~isnan(cent_cells), 2) & (DataInfo(:, 1) == 0 | DataInfo(:, 2) == 0), 2),...
    %     25,  colors(i, :), 'filled');
end
plot([0 800], 300*ones(1, 2), '--k');
plot(400*ones(1, 2), [0 600], '--k');
% xticks(0:3);
% xticklabels({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});
xlabel('Width (pixel)');
ylabel('Height (pixel)');
xlim([300 500]);
ylim([150 400]);
title('OFF cell center by retina')
%% Quatify
retina_cent_count = cellfun(@(x) size(x, 1), retina_cent);
retina_cent(retina_cent_count<2) = [];
all_retina_cent = [];
paired_retina_cent = [];
retina_cent_count = cellfun(@(x) size(x, 1), retina_cent);
for i = 1:length(retina_cent)
    ccent = retina_cent{i};
    all_retina_cent = [all_retina_cent; ccent];
    pids = nchoosek(1:retina_cent_count(i), 2);
    for j = 1:size(pids, 1)
        paired_retina_cent = [paired_retina_cent; ccent(pids(j, 1), :), ccent(pids(j, 2), :)];
    end
end
%%
num_sample = 10000;
paired_distance = nan(num_sample, 2);
pixel_to_um = 2.5;
rng('shuffle');
for i = 1:num_sample
    sids = randsample(1:size(all_retina_cent, 1), 2, false);
    paired_distance(i, 2) = sqrt(sum((all_retina_cent(sids(1), :)-all_retina_cent(sids(2), :)).^2))*pixel_to_um;
    sids = randsample(1:size(paired_retina_cent, 1), 1);
    paired_distance(i, 1) = sqrt(sum((paired_retina_cent(sids, 1:2)-paired_retina_cent(sids, 3:4)).^2))*pixel_to_um;
end
paired_distance = paired_distance + 0.1*randn(size(paired_distance));
bins = linspace(min(paired_distance(:)), max(paired_distance(:)), 30);
figure; hold on
[f, x, fl, fu] = ecdf(paired_distance(:, 1));
shadePlot(x, f, fu-f, 'k', 0.5);
[f, x, fl, fu] = ecdf(paired_distance(:, 2));
shadePlot(x, f, fu-f, 'r', 0.5);
xlabel('Paired distance (um)');
ylabel('Probability');
ylim([0 1])


x = [1*ones(length(unique(paired_distance(:, 1))), 1); 2*ones(length(unique(paired_distance(:, 2))), 1)];
y = [unique(paired_distance(:, 1)); unique(paired_distance(:, 2))];
figure; hold on
swarmchart(x, y, 25, 'k', 'filled');
plot([0.7 1.3], median(unique(paired_distance(:, 1)))*ones(1, 2), 'k');
plot([1.7 2.3], median(unique(paired_distance(:, 2)))*ones(1, 2), 'k');
xticks(1:2)
xticklabels({'within', 'between'});
xlabel('Paired condition');
ylabel('Paired distance (um)');

figure; hold on
h1 = histogram(paired_distance(:, 1), bins);
h1.Normalization = 'probability';
h1.FaceColor = 'k';
h2 = histogram(paired_distance(:, 2), bins);
h2.Normalization = 'probability';
h2.FaceColor = 'r';
h2.FaceAlpha = 0.5;

%%
diame_ids = loadedFile.dim2_diameters;
speed_ids = loadedFile.dim3_speeds;
direc_ids = loadedFile.dim1_moving_direction;

%% Make sure temporal filter is matching their cell type
cids = DataInfo(:, 1)*2+DataInfo(:, 2);
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
t = WinT(1):1/Fz:WinT(end);
type_names = {'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'};
avg_tf = nan(length(ucids), size(Data_TF, 2));
figure; hold on
for i = ucids(:)'
    switch i
        case 2 % ON-Nasal :10
            line_color = [92 204 206]/255;
            subfig = 3;
        case 3 % ON-Temporal :11
            line_color = [255 139 139]/255;
            subfig = 4;
        case 0 % OFF-Nasal :0
            line_color = [60 91 89]/255;
            subfig = 1;
        case 1 % OFF-Temporal :1
            line_color = [136 74 75]/255;
            subfig = 2;
    end
    ca = Data_TF(cids == i, :);
    ca = ca./max(abs(ca), [], 2);
    avg_tf(i+1, :) = mean(ca, 1);
    subplot(2, 2, subfig); hold on
    plot(repmat(t(2:end)', 1, size(ca, 1)), ca', 'Color', line_color);
    plot(t(2:end), mean(ca, 1), 'k', 'LineWidth',2)
    xlim([t(2), t(end)])
    xlabel('Time (s)');
    ylabel('Effective contrast');
    title(type_names{subfig});
end

%%

figure;
subplot(1, 2, 1);
ca = Data_TF(DataInfo(:, 1) == 1, :);
ca = ca./max(abs(ca), [], 2);
plot(t(2:end), ca);
hold on
plot(t(2:end), mean(ca, 1), 'k', 'LineWidth',2)
xlabel('Time (s)');
ylabel('STA contrast');
ylim([-1 1]);
title('ON temporal kernel')
box off
subplot(1, 2, 2)
ca = Data_TF(DataInfo(:, 1) == 0, :);
ca = ca./max(abs(ca), [], 2);
plot(t(2:end), ca);
hold on
plot(t(2:end), mean(ca, 1), 'k', 'LineWidth',2)
xlabel('Time (s)');
ylabel('STA contrast');
ylim([-1 1]);
box off
title('OFF temporal kernel')
%%
% save('./Results/SummaryData_041524.mat');
% save('./Results/SummaryData_042324.mat'); % optimal parameters
% save('./Results/SummaryData_051724.mat'); 
%% Middle point loading
% load('./Results/SummaryData_041524.mat');
% save('./Results/SummaryData_051724.mat');
load('./Results/SummaryData_051724.mat'); 
%%
keyboard;
%% plot the Response of moving spot and examine cell type bases
close all
spot_ids = ~cellfun(@isempty, Data_Spot);
spot_ids = find(spot_ids & ~isnan(DataInfo(:, 1)'));
size_count = []; % use for set the minimum length
t_length = 369;
std_thr = 2.58;
baseline_time = 0.3;
DisplaySquare = 4;
summary_info = [];
summary_trace = [];
for i = 1:length(direc_ids)
    figure;
    for r = 1:length(spot_ids)
        spot_id = spot_ids(r);
        Data = Data_Spot{spot_id};
        if isempty(Data)
            continue
        end
        size_count = [size_count; i, r, size(Data)];
        csignal = reshape(medfilt1(mean(Data(i, :, 1, :, 1:round(baseline_time*Fz)), 4), 5, [], 5), [], 1);
        signal_mean = mean(csignal);
        signal_std = std(csignal);
        type_id = 10*DataInfo(spot_id, 1)+DataInfo(spot_id, 2);
        switch type_id
            case 10 % ON-Nasal
                line_color = [92 204 206]/255;
            case 11 % ON-Temporal
                line_color = [255 139 139]/255;
            case 0 % OFF-Nasal
                line_color = [60 91 89]/255;
            case 1 % OFF-Temporal
                line_color = [136 74 75]/255;
        end
        t = (0:(size(Data, 5)-1))/Fz;
        ii = 1;
        for k = length(diame_ids)-DisplaySquare:length(diame_ids)
            for j = length(speed_ids)-DisplaySquare:length(speed_ids)
                subplot(DisplaySquare+1, DisplaySquare+1, ii); hold on
                condition_signal =squeeze(Data(i, k, j, :, :));
                cal_type = 1;
                switch cal_type
                    case 1
                        avg_signal = medfilt1(mean(condition_signal, 1), 5);
                    case 2
                        avg_signal = std(medfilt1(condition_signal, 5, [], 2), [], 1);
                    case 3
                        avg_signal = medfilt1(mean(condition_signal, 1), 5);
                        avg_signal = cumsum(diff([avg_signal(1) avg_signal]));
                end
                % plot(t, condition_signal, 'Color', 0.5*ones(1, 3));
                summary_info = [summary_info; i, r, type_id, k, j, DataInfo(spot_id, 3), spot_id];
                summary_trace = [summary_trace; avg_signal(1:t_length)];
                plot(t, avg_signal, 'Color', line_color);
                cids = find(abs(avg_signal-signal_mean) > signal_std*std_thr);

                stop_time = mean(Data_MovStop{spot_id}(i, j, j, :)) ;
                if ~isempty(cids)
                    plot(t(cids(1)), avg_signal(cids(1)), 'xr');
                    resp_onset(i, k, j) = stop_time - cids(1);
                end
                plot(t(round(stop_time))*ones(1, 2), [0 100], 'b');
                % xlim([0, 0.8]);
                if k == 3
                    title(sprintf('Speeds:%d(um/s)', speed_ids(j)))
                    Label_spds{j} = sprintf('%d (um/s)', speed_ids(j));
                end
                if j == 1
                    if k ==1
                        ylabel(sprintf('Diameters:%d (um)', diame_ids(k)));
                    else
                        ylabel(sprintf('%d (um)', diame_ids(k)));
                    end
                    Label_dias{k} = sprintf('%d (um)', diame_ids(k));
                elseif j == 2
                    ylabel('Firing rate (spike/s)')
                end
                % ylim([0 150]);
                ii = ii + 1;
            end
        end
        sgtitle(sprintf('%s Directions:%d', direc_ids(i)));
    end
end

%% Get template traces for each type in each condition under moving spots
close all
t = (0:(size(Data, 5)-1))/Fz;
cell_type_list = [0 1 10 11];
diame_exam_ids = 6:7;
speed_exam_ids = 3:5;
template_traces = nan(length(direc_ids)*length(cell_type_list)*length(diame_exam_ids)*length(speed_exam_ids), t_length);
template_info = nan(length(direc_ids)*length(cell_type_list)*length(diame_exam_ids)*length(speed_exam_ids), 4);
jj = 1;
for i = 1:length(direc_ids)
    figure;
    ii = 1;
    for k = diame_exam_ids
        for j = speed_exam_ids
            subplot(length(diame_exam_ids), length(speed_exam_ids), ii); hold on
            for q = 1:length(cell_type_list)
                cids = summary_info(:, 1) == i & summary_info(:, 4) == k &...
                    summary_info(:, 5) == j & summary_info(:, 3) == cell_type_list(q);
                avg_signal = mean(summary_trace(cids, :), 1);
                std_signal = std(summary_trace(cids, :), [], 1)/sqrt(sum(cids));
                switch cell_type_list(q)
                    case 10 % ON-Nasal
                        line_color = [92 204 206]/255;
                    case 11 % ON-Temporal
                        line_color = [255 139 139]/255;
                    case 0 % OFF-Nasal
                        line_color = [60 91 89]/255;
                    case 1 % OFF-Temporal
                        line_color = [136 74 75]/255;
                end
                shadePlot(t(1:t_length), avg_signal, std_signal*1.96, line_color);
                %plot(t(1:t_length), avg_signal, 'Color', line_color);
                template_traces(jj, :) = avg_signal;
                template_info(jj, :) = [i, q, k, j];
                jj = jj +1;
            end

            stop_time = mean(Data_MovStop{spot_id}(i, j, j, :)) ;
            plot(t(round(stop_time))*ones(1, 2), [0 100], '--k');

            xlim([0, t(sum(~isnan(avg_signal)))]);
            ylim([0 150]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)');

            if k == 3
                title(sprintf('Speeds:%d(um/s)', speed_ids(j)))
                Label_spds{j} = sprintf('%d (um/s)', speed_ids(j));
            end
            if j == 1
                if k ==1
                    ylabel(sprintf('Diameters:%d (um)', diame_ids(k)));
                else
                    ylabel(sprintf('%d (um)', diame_ids(k)));
                end
                Label_dias{k} = sprintf('%d (um)', diame_ids(k));
            elseif j == 2
                ylabel('Firing rate (spike/s)')
            end
            % ylim([0 150]);
            ii = ii + 1;
        end
    end
    legend({'OFFnasal', '', 'OFFtemporal', '', 'ONnasal', '', 'ONtemporal', ''});
    sgtitle(sprintf('Directions:%d',  direc_ids(i)));
end


%% Calculate the orientation and see if the results matching in correlation
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
cids = DataInfo(:, 3); % retina specific
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
colors = lines(length(ucids));
DisplaySquare = 1; % the last DisplaySquare + 1
contrast_type = [0 1; 10 11];% [0 1; 10 11]
up_sample_fac = 10;
IsDisplay = 0;
pixel_size = 2.5;
diame_exam_ids = 6:7;
speed_exam_ids = 3:5;
close all
for i = 1:length(direc_ids)
    ii = 1;
    posi_info = [];
    delay_info = [];
    type_info = [];
    for u = ucids(:)'
        clc
        fprintf('progress... %d/%d \n', u, ucids(end));
        for k =  diame_exam_ids % length(diame_ids)-DisplaySquare:length(diame_ids)
            for j = speed_exam_ids % length(speed_ids)-DisplaySquare:length(speed_ids)
                for q = 2
                    cids = find(summary_info(:, 1) == i & summary_info(:, 4) == k &...
                        summary_info(:, 5) == j & ismember(summary_info(:, 3), contrast_type(q, :)) &...
                        summary_info(:, 6) == u);
                    cell_ids = summary_info(cids, 2);
                    ncell = length(cids);
                    if ncell < 2
                        continue
                    end
                    pair_ids = nchoosek(cids, 2);
                    npair = size(pair_ids, 1);
                    for p = 1:npair

                        s1 = summary_trace(pair_ids(p, 1), :);
                        s2 = summary_trace(pair_ids(p, 2), :);
                        exist_time = all(~isnan([s1; s2]), 1);
                        s1 = s1(exist_time);
                        s2 = s2(exist_time);
                        t = (0:length(s1)-1)/Fz;

                        s1_up = SimpleSignalUpSample(s1(exist_time), up_sample_fac);
                        s2_up = SimpleSignalUpSample(s2(exist_time), up_sample_fac);
                        t_up = (0:length(s1_up)-1)/(Fz*up_sample_fac);

                        [c, lags] = xcorr(s1_up(:), s2_up(:));
                        [~, iLag] = max(c);
                        s3_up = circshift(s2_up, [0 iLag]);
                        shift_time = lags(iLag)/(Fz*up_sample_fac);


                        pos_1 = cent_cells(summary_info(pair_ids(p, 1), 7), :);
                        pos_2 = cent_cells(summary_info(pair_ids(p, 2), 7), :);
                        posi_info = [posi_info; pos_2, pos_1];
                        type_info = [type_info; summary_info(pair_ids(p, 1), 2:3)];
                        delay_info = [delay_info; shift_time, diame_ids(k), speed_ids(j), u, k, j, p];

                        if IsDisplay
                            close all
                            figure; hold on
                            plot(t, s1, '--k');
                            plot(t, s2, '--m');
                            plot(t_up, s1_up, 'k');
                            plot(t_up, s2_up, 'm');
                            plot(t_up, s3_up, 'b');
                            title(sprintf('Lag: %.03G (s)', shift_time));
                            keyboard;
                        end
                    end

                end
            end
        end
    end


    % Example the delay and time shift
    Colors = [60 91 89;
        136 74 75;
        92 204 206;
        255 139 139]/255;

    figure;
    subplot(1, 2, 1); hold on
    x = (posi_info(:, 1)-posi_info(:, 3))*pixel_size;
    y = delay_info(:, 1).*delay_info(:, 3);
    [yFit, beta] = SimpleRegression2D(x, y);
    [~, sids] = sort(x);
    gscatter(x, y, type_info(:, 2));
    plot(x(sids), yFit(sids), '--r');
    xlabel('Position distance (um)');
    ylabel('Moving spot delay distance (um)');
    % plot([-150 150], [150 -150], '--k');
    plot([-150 150], [-150 150], '--k');
    title(sprintf('Width difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));

    subplot(1, 2, 2); hold on
    x = (posi_info(:, 2)-posi_info(:, 4))*pixel_size;
    y = delay_info(:, 1).*delay_info(:, 3);
    [yFit, beta] = SimpleRegression2D(x, y);
    [~, sids] = sort(x);
    % scatter(x, y, 25, 'k', 'filled');
    gscatter(x, y, type_info(:, 2));
    plot(x(sids), yFit(sids), '--r');
    plot([-150 150], [-150 150], '--k');
    xlabel('Position distance (um)');
    ylabel('Moving spot delay distance (um)');
    title(sprintf('Height difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));

    sgtitle(sprintf('Direction %d', direc_ids(i)));

    % subplot(1, 4, 3); hold on
    % x = sqrt(sum((posi_info(:, 1:2)-posi_info(:, 3:4)).^2, 2))*pixel_size;
    % y = delay_info(:, 1).*delay_info(:, 3);
    % [yFit, beta] = SimpleRegression2D(x, y);
    % [~, sids] = sort(x);
    % scatter(x, y, 25, 'k', 'filled');
    % plot(x(sids), yFit(sids), '--r');
    % xlabel('Position distance (um)');
    % ylabel('Moving spot delay distance (um)');
    % title(sprintf('Spatial difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));
    %
    % subplot(1, 4, 4); hold on
    % x = (sqrt(sum(([800 600] -posi_info(:, 3:4)).^2, 2)))*pixel_size;
    % y = delay_info(:, 1).*delay_info(:, 3);
    % [yFit, beta] = SimpleRegression2D(x, y);
    % [~, sids] = sort(x);
    % scatter(x, y, 25, 'k', 'filled');
    % plot(x(sids), yFit(sids), '--r');
    % xlabel('Position distance (um)');
    % ylabel('Moving spot delay distance (um)');
    %
    % title(sprintf('To-Center difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));

end

%% [Section C] Alignment to the global template
cids = DataInfo(:, 3); % retina specific
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
colors = lines(length(ucids));
DisplaySquare = 1; % the last DisplaySquare + 1
up_sample_fac = 10;
IsDisplay = 0;
pixel_size = 2.5;
diame_exam_ids = 6:7;
speed_exam_ids = 3:5;
close all
di = 2 % 1:length(direc_ids)
ii = 1;
posi_info = [];
delay_info = [];

figure;
for u = ucids(:)'
    clc
    fprintf('progress... %d/%d \n', u, ucids(end));
    for k =  1:length(diame_exam_ids) % length(diame_ids)-DisplaySquare:length(diame_ids)
        for j = 1:length(speed_exam_ids) % length(speed_ids)-DisplaySquare:length(speed_ids)
            subplot()
            for q = [1:2]% [1:2] 1:length(cell_type_list)
                cids = find(summary_info(:, 1) == di & summary_info(:, 4) == diame_exam_ids(k) &...
                    summary_info(:, 5) == speed_exam_ids(j) & ismember(summary_info(:, 3), cell_type_list(q)) &...
                    summary_info(:, 6) == u);
                ncell = length(cids);

                tid = template_info(:, 1) == di & template_info(:, 2) == q &...
                    template_info(:, 3) == diame_exam_ids(k) & template_info(:, 4) == speed_exam_ids(j);
                for p = 1:ncell
                    s1 = template_traces(tid,  :);
                    s2 = summary_trace(cids(p), :);
                    exist_time = all(~isnan([s1; s2]), 1);
                    s1 = s1(exist_time);
                    s2 = s2(exist_time);
                    t = (0:length(s1)-1)/Fz;

                    s1_up = SimpleSignalUpSample(s1(exist_time), up_sample_fac);
                    s2_up = SimpleSignalUpSample(s2(exist_time), up_sample_fac);
                    t_up = (0:length(s1_up)-1)/(Fz*up_sample_fac);

                    [c, lags] = xcorr(s1_up(:), s2_up(:));
                    [~, iLag] = max(c);
                    s3_up = circshift(s2_up, [0 iLag]);
                    shift_time = lags(iLag)/(Fz*up_sample_fac);

                    pos_1 = [400 300];
                    pos_2 = cent_cells(summary_info(cids(p), 7), :);
                    posi_info = [posi_info; pos_2, pos_1];
                    delay_info = [delay_info; shift_time, diame_ids(k), speed_ids(j),...
                        u, diame_exam_ids(k), speed_exam_ids(j), p, summary_info(cids(p), 7)];

                    if IsDisplay
                        close all
                        figure; hold on
                        plot(t, s1, '--k');
                        plot(t, s2, '--m');
                        plot(t_up, s1_up, 'k');
                        plot(t_up, s2_up, 'm');
                        plot(t_up, s3_up, 'b');
                        title(sprintf('Lag: %.03G (s)', shift_time));
                        keyboard;
                    end
                end

            end
        end
    end
end


figure;
subplot(1, 2, 1); hold on
x = (posi_info(:, 1)-posi_info(:, 3))*pixel_size;
y = delay_info(:, 1).*delay_info(:, 3);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
if ismember(direc_ids(di), [0 90])
    plot([-150 150], [150 -150], '--k');
else
    plot([-150 150], [-150 150], '--k');
end
xlabel('Position distance (um)');
ylabel('Moving spot delay distance (um)');
title(sprintf('Width difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));

subplot(1, 2, 2); hold on
x = (posi_info(:, 2)-posi_info(:, 4))*pixel_size;
y = delay_info(:, 1).*delay_info(:, 3);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
plot([-150 150], [-150 150], '--k');
xlabel('Position distance (um)');
ylabel('Moving spot delay distance (um)');
title(sprintf('Height difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));

% subplot(1, 4, 3); hold on
% x = sqrt(sum((posi_info(:, 1:2)-posi_info(:, 3:4)).^2, 2))*pixel_size;
% y = delay_info(:, 1).*delay_info(:, 3);
% [yFit, beta] = SimpleRegression2D(x, y);
% [~, sids] = sort(x);
% scatter(x, y, 25, 'k', 'filled');
% plot(x(sids), yFit(sids), '--r');
% xlabel('Position distance (um)');
% ylabel('Moving spot delay distance (um)');
% title(sprintf('Spatial difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));
%
% subplot(1, 4, 4); hold on
% x = (sqrt(sum(([800 600] -posi_info(:, 1:2)).^2, 2)))*pixel_size;
% y = delay_info(:, 1).*delay_info(:, 3);
% [yFit, beta] = SimpleRegression2D(x, y);
% [~, sids] = sort(x);
% scatter(x, y, 25, 'k', 'filled');
% plot(x(sids), yFit(sids), '--r');
% xlabel('Position distance (um)');
% ylabel('Moving spot delay distance (um)');
% title(sprintf('To-Center difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));

sgtitle(sprintf('Direction %d', direc_ids(di)));

%% Visualize (corrected shift) responses of ON cells according to the position of OFF cells
% (1) use 1 to 1 ratio (for the correction according to the position
% (2) shift of ON reponses to spot will be based on the center of OFF cell
%     in the same retina
% (3) (control) Plot OFF response correlation, and it should have little
%     effect
% (4) expecting ON nasal and ON temporal will be different due to the
%     in-center contrast (heterogenity)
% (5) Instead of comparing nasal and temporal, we can also see the relation
%     betweeen gaussian-likenss and correct off-distance

cids = DataInfo(:, 3); % retina specific
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
colors = lines(length(ucids));
up_sample_fac = 10;
IsDisplay =0;
pixel_size = 2.5;
diame_exam_ids = 6:7;
speed_exam_ids = 3:5;
anchor_pos = [400 300];
max_cell_per_type = 5;
cell_type_list = [0, 1, 10, 11];
is_corrected = 1;
% close all
for i =  1:length(direc_ids)
    ii = 1;
    switch direc_ids(i)
        case 0
            pos_to_shiftlag = 1;
            shift_pos_ind = 1;
        case 90
            pos_to_shiftlag = -1;
            shift_pos_ind = 2;
        case 180
            pos_to_shiftlag = -1;
            shift_pos_ind = 1;
    end
    gather_responses = nan(max_cell_per_type*length(diame_exam_ids)*length(speed_exam_ids)*4, size(summary_trace, 2)*up_sample_fac);
    gather_info = nan(max_cell_per_type*length(diame_exam_ids)*length(speed_exam_ids)*4, 3);
    gi = 1;
    for u = ucids(:)'
        clc
        fprintf('progress... %d/%d \n', u, ucids(end));
        for k =  1:length(diame_exam_ids) % length(diame_ids)-DisplaySquare:length(diame_ids)
            for j = 1:length(speed_exam_ids) % length(speed_ids)-DisplaySquare:length(speed_ids)
                % get OFF cells in the recording, and align ON cells to
                % cell_type_list = [0 1 10 11] for summary_info(:, 3)

                cids = find(summary_info(:, 1) == i & summary_info(:, 4) == diame_exam_ids(k) &...
                    summary_info(:, 5) == speed_exam_ids(j) & summary_info(:, 6) == u);
                type_ids = summary_info(cids, 3);
                if ~any(ismember(type_ids, [0 1])) % skip if there is no off cells
                    continue
                end
                % find the OFF center, and all cells will be shifted
                off_ids = cids(ismember(type_ids, [0 1]));
                off_center_pos = mean(cent_cells(summary_info(off_ids, 7), :), 1);
                ncell = length(cids);
                for p = 1:ncell
                    s2 = summary_trace(cids(p), :);
                    exist_time = all(~isnan(s2), 1);
                    s2 = s2(exist_time);
                    t = (0:length(s2)-1)/Fz;
                    s2_up = SimpleSignalUpSample(s2, up_sample_fac);
                    t_up = (0:length(s2_up)-1)/(Fz*up_sample_fac);
                    rf_offset = sqrt(sum((off_center_pos(shift_pos_ind) - anchor_pos(shift_pos_ind)).^2));
                    offset_time = rf_offset*pixel_size/speed_ids(j);
                    iLag = round(pos_to_shiftlag*Fz*up_sample_fac*offset_time + length(s2_up));
                    s2_up_shift = circshift(s2_up, [0 iLag]);

                    if is_corrected
                        gather_responses(gi, 1:length(s2_up)) = s2_up_shift;
                    else
                        gather_responses(gi, 1:length(s2_up)) = s2_up;
                    end
                    gather_info(gi, :) = [type_ids(p), diame_exam_ids(k), speed_exam_ids(j)];
                    gi = gi +1;
                    %plot(t_up, s2_up, 'Color', line_color);
                    % plot(t_up, s2_up_shift, 'Color', line_color);
                    if IsDisplay
                        close all
                        figure; hold on
                        plot(t_up, s2_up, 'k');
                        plot(t_up, s2_up_shift, 'm');
                        title(sprintf('displacement: %3G (um)', rf_offset));
                        keyboard;
                    end
                end
            end
        end
    end


    figure;
    ii = 1;
    for k =  1:length(diame_exam_ids) % length(diame_ids)-DisplaySquare:length(diame_ids)
        for j = 1:length(speed_exam_ids) % length(speed_ids)-DisplaySquare:length(speed_ids)
            subplot(length(diame_exam_ids), length(speed_exam_ids), ii); hold on
            ii = ii + 1;
            for q = 1:length(cell_type_list)
                switch cell_type_list(q)
                    case 10 % ON-Nasal
                        line_color = [92 204 206]/255;
                    case 11 % ON-Temporal
                        line_color = [255 139 139]/255;
                    case 0 % OFF-Nasal
                        line_color = [60 91 89]/255;
                    case 1 % OFF-Temporal
                        line_color = [136 74 75]/255;
                end
                cids = gather_info(:, 1) == cell_type_list(q) & gather_info(:, 2) == diame_exam_ids(k) &...
                    gather_info(:, 3) == speed_exam_ids(j);
                avg = mean(gather_responses(cids, :), 1);
                sem = std(gather_responses(cids, :), [], 1)/sqrt(sum(cids));
                exist_time = all(~isnan([avg; sem]), 1);
                avg = avg(exist_time);
                sem = sem(exist_time);
                t = (0:length(avg)-1)/(Fz*up_sample_fac);
                shadePlot(t, avg, sem*1.96, line_color);
            end

            xlim([0, t(sum(~isnan(avg)))]);
            ylim([0 150]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)');

            if k == 1
                title(sprintf('Speeds:%d(um/s)', speed_ids(speed_exam_ids(j))))
                Label_spds{j} = sprintf('%d (um/s)', speed_ids(speed_exam_ids(j)));
            end
            if j == 1
                if k ==1
                    ylabel(sprintf('Diameters:%d (um)', diame_ids(diame_exam_ids(k))));
                else
                    ylabel(sprintf('%d (um)', diame_ids(diame_exam_ids(k))));
                end
                Label_dias{k} = sprintf('%d (um)', diame_ids(diame_exam_ids(k)));
            elseif j == 2
                ylabel('Firing rate (spike/s)')
            end
        end
    end
    if is_corrected
        sgtitle(sprintf('Direction:%d (corrected)', direc_ids(i)));
    else
        sgtitle(sprintf('Direction:%d (uncorrected)', direc_ids(i)));
    end
end

%% Shift temporal filter based on responses from white noise
% *** run Section C first for a particular direction
% check if the temporal filter stretching can recapitulate the moving spot
% responses shift
% TF shift
% (1) get avg_tf (from previous code for getting temporal filter)
t = WinT(1):1/Fz:WinT(end);
up_sample_fac = 10;
IsDisplay = 0;
cell_type_list = [0, 1, 10, 11];
tf_delay_time = nan(size(Data_TF, 1), 2);
for i = 1:size(Data_TF, 1)
    type_id = 10*DataInfo(i, 1)+DataInfo(i, 2);
    s1 = avg_tf(cell_type_list==type_id, :);
    s2 = Data_TF(i, :);
    s2 = s2./max(abs(s2), [], 2);

    exist_time = all(~isnan([s1; s2]), 1);
    s1 = s1(exist_time);
    s2 = s2(exist_time);

    s1_up = SimpleSignalUpSample(s1(exist_time), up_sample_fac);
    s2_up = SimpleSignalUpSample(s2(exist_time), up_sample_fac);
    t_up = SimpleSignalUpSample(t(2:end), up_sample_fac);

    [c, lags] = xcorr(s1_up(:), s2_up(:));
    [~, iLag] = max(c);
    s3_up = circshift(s2_up, [0 iLag]);
    shift_time = lags(iLag)/(Fz*up_sample_fac);
    tf_delay_time(i, :) = [shift_time, type_id];
    if IsDisplay
        close all
        figure; hold on
        plot(t(2:end), s1, '--k');
        plot(t(2:end), s2, '--m');
        plot(t_up, s1_up, 'k');
        plot(t_up, s2_up, 'm');
        plot(t_up, s3_up, 'b');
        title(sprintf('Lag: %.03G (s)', shift_time));
        keyboard;
    end
end

%
figure;
subplot(1, 2, 1); hold on
x = -tf_delay_time(delay_info(:, 8), 1);
y = delay_info(:, 1).*delay_info(:, 3);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('TF delay (s)');
ylabel('Moving spot delay distance (um)');
title(sprintf('Position difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:)))); % Width(1, 3) Height (2)

subplot(1, 2, 2); hold on
x = -tf_delay_time(delay_info(:, 8), 1);
y = delay_info(:, 1);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('TF delay (s)');
ylabel('Moving spot delay (s)');
title(sprintf('Position difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));


sgtitle(sprintf('Direction %d', direc_ids(di)));

%% Stretching the responses based on temporal filter from white noise
% *** run Section C first for a particular direction
alignment_type = 'xcorr'; % 'peak' or 'xcorr'
t = WinT(1):1/Fz:WinT(end);
up_sample_fac = 10;
IsDisplay = 0;
cell_type_list = [0, 1, 10, 11];
Stretching_info = nan(size(Data_TF, 1), 6);
%tf_delay_time = nan(size(Data_TF, 1), 2);
for i = 1:size(Data_TF, 1)
    type_id = 10*DataInfo(i, 1)+DataInfo(i, 2);
    s1 = avg_tf(cell_type_list==type_id, :);
    if ismember(type_id, [0 1])
        s4 = mean(avg_tf(ismember(cell_type_list, [0 1]), :), 1);
    else
        s4 = mean(avg_tf(ismember(cell_type_list, [10 11]), :), 1);
    end
    s2 = Data_TF(i, :);
    s2 = s2./max(abs(s2), [], 2);

    exist_time = all(~isnan([s1; s2]), 1);
    s1 = s1(exist_time);
    s2 = s2(exist_time);
    s4 = s4(exist_time);

    s1_up = SimpleSignalUpSample(s1(exist_time), up_sample_fac);
    s2_up = SimpleSignalUpSample(s2(exist_time), up_sample_fac);
    s4_up = SimpleSignalUpSample(s4(exist_time), up_sample_fac);
    t_up = SimpleSignalUpSample(t(2:end), up_sample_fac);
    % first align them to the peak and stretch at the peak position
    switch lower(alignment_type)
        case 'peak'
            if ismember(type_id, [0 1])
                [~, s1_c] = max(-s1_up);
                [~, s2_c] = max(-s2_up);
            else
                [~, s1_c] = max(s1_up);
                [~, s2_c] = max(s2_up);
            end
            iLag = length(s1_up)+(s1_c-s2_c);

        case 'xcorr'
            [c, lags] = xcorr(s1_up(:), s2_up(:));
            [~, iLag] = max(c);
            [c, lags_3] = xcorr(s4_up(:), s2_up(:));
            [~, iLag_3] = max(c);
    end
    s3_up = circshift(s2_up, [0 iLag]);
    s5_up = circshift(s2_up, [0 iLag_3]);
    if iLag > length(s2_up)
        s3_up(1:(iLag-length(s2_up))) = s2_up(1);
        s5_up(1:(iLag_3-length(s2_up))) = s2_up(1);
    else
        s3_up(end-(length(s2_up)-iLag):end) = s2_up(end);
        s5_up(end-(length(s2_up)-iLag_3):end) = s2_up(end);
    end
    shift_time = lags(iLag)/(Fz*up_sample_fac);
    shift_time_3 = lags_3(iLag_3)/(Fz*up_sample_fac);
    %tf_delay_time(i, :) = [shift_time, type_id];
    if IsDisplay
        close all
        figure;
        subplot(1, 2, 1); hold on
        plot(t(2:end), s1, '--k');
        plot(t(2:end), s2, '--m');
        plot(t_up, s1_up/max(abs(s1_up)), 'k');
        plot(t_up, s2_up, 'm');
        plot(t_up, s3_up, 'b');
        legend({'Template', 'Source', 'Norm. template', 'Up sampled source', 'Shifted'});
        title(sprintf('Lag: %.03G (s)', shift_time));
    end

    % Perform stretching to the aligned
    [~, max_id] = max(s3_up.*s1_up);
    t_stretching = t_up-t_up(max_id);
    t_axis_st = stretchKernel(t_stretching, s3_up, 1.3);
    t_axis_st = t_axis_st + t_up(max_id);

    [~, max_id] = max(s5_up.*s1_up);
    t_stretching_3 = t_up-t_up(max_id);


    % figure;  hold on
    % plot(t_up, s3_up, 'k');
    % plot(t_axis_st, s3_up, 'b');

    % Optimization to find the best stretch factor
    lb = 0.5;
    ub = 2;

    options = optimset('Display', 'off', 'TolX', 1e-6);
    best_stretch = fminbnd(@(x) stretchingcostFunction(x, s3_up, s1_up, t_stretching), lb, ub, options);
    best_stretch_3 = fminbnd(@(x) stretchingcostFunction(x, s5_up, s4_up, t_stretching_3), lb, ub, options);

    % interpolate to the same axis points (use the smallest coverage one)
    r = -stretchingcostFunction(best_stretch, s3_up, s1_up, t_stretching);

    t_axis_st = stretchKernel(t_stretching, s3_up, best_stretch);
    t_axis_st = t_axis_st + t_up(max_id);
    Stretching_info(i, :) = [best_stretch, r, type_id, best_stretch_3, shift_time, shift_time_3];

    if IsDisplay
        subplot(1, 2, 2); hold on
        plot(t_up, s1_up/max(abs(s1_up)), 'k');
        plot(t_up, s3_up, '--b');
        plot(t_axis_st, s3_up, 'b');
        legend({'Template', 'Source', 'Stretched'});
        title(sprintf('stretching factor: %.02G (r=%.03G)', best_stretch, r));
        keyboard;
    end
end

%
figure;
subplot(1, 2, 1); hold on
x = Stretching_info(delay_info(:, 8), 1);
y = delay_info(:, 1).*delay_info(:, 3);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('Stretching factor');
ylabel('Moving spot delay distance (um)');
title(sprintf('Position difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:)))); % Width(1, 3) Height (2)

subplot(1, 2, 2); hold on
x = Stretching_info(delay_info(:, 8), 1);
y = delay_info(:, 1);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('Stretching factor');
ylabel('Moving spot delay (s)');
title(sprintf('Position difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:))));

sgtitle(sprintf('Direction %d', direc_ids(di)));
%%
figure; hold on
x = Stretching_info(:, 4);
y = Stretching_info(:, 6);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('Stretching factor');
ylabel('TF delay (s)');
title(sprintf('Correlation betweeen kernel parameters (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:)))); % Width(1, 3) Height (2)
%% Summary plot TF kernel analysis with shift and stretching of % EXAMINE by RETINA (BATCH EFFECT) unfinished
figure;
% subplot(1, 4, 1); hold on
x = nan(size(Stretching_info, 1), 1);
y = nan(size(Stretching_info, 1), 1);
for i = 1:size(Stretching_info)
    x(i) = find(cell_type_list == Stretching_info(i, 3));
    y(i) = Stretching_info(i, 4);
end
swarmchart(x, y, 25, 'k', 'filled');
linelength = 0.6;
for i = 1:4
    plot([i-linelength/2 i+linelength/2],...
        median(Stretching_info(Stretching_info(:, 3) == cell_type_list(i), 4))*ones(1, 2), 'k');
end
xlabel('cell types');
ylabel('Stretching factor');
xticks(1:4)
xticklabels({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});

subplot(1, 3, 1); hold on
x = nan(size(Stretching_info, 1), 1);
y = nan(size(Stretching_info, 1), 1);
for i = 1:size(Stretching_info)
    x(i) = find(cell_type_list == Stretching_info(i, 3));
    y(i) = Stretching_info(i, 6);
end
swarmchart(x, y, 25, 'k', 'filled');
for i = 1:4
    plot([i-linelength/2 i+linelength/2],...
        median(Stretching_info(Stretching_info(:, 3) == cell_type_list(i), 6))*ones(1, 2), 'k');
end
xlabel('cell types');
ylabel('Shift time (s)');
xticks(1:5)
xticklabels({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});

subplot(1, 3, 2); hold on
x = nan(size(Stretching_info, 1), 1);
y = nan(size(Stretching_info, 1), 1);
for i = 1:size(Stretching_info)
    x(i) = find(cell_type_list == Stretching_info(i, 3));
    y(i) = Stretching_info(i, 6);
end
gscatter(x+linelength*(rand(size(x, 1), 1)-0.5), y, DataInfo(:, 3));
for i = 1:4
    plot([i-linelength/2 i+linelength/2],...
        median(Stretching_info(Stretching_info(:, 3) == cell_type_list(i), 6))*ones(1, 2), 'k');
end
xlabel('cell types');
ylabel('Shift time (s)');
xticks(1:5)
xticklabels({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});

subplot(1, 3, 3); hold on
corrected_y = nan(size(Stretching_info, 1), 1);
x = nan(size(Stretching_info, 1), 1);
for i = 1:size(Stretching_info)
    cids = find(DataInfo(:, 3) == DataInfo(i, 3));
    sids = DataInfo(cids, 1) == 0;
    if sum(sids) == 0
        continue
    else
        x(i) = find(cell_type_list == Stretching_info(i, 3));
        corrected_y(i) = Stretching_info(i, 6) - mean(Stretching_info(cids(sids), 6));
    end
end
gids = DataInfo(~isnan(corrected_y), 3);
x = x(~isnan(corrected_y));
corrected_y = corrected_y(~isnan(corrected_y));
gscatter(x+linelength*(rand(size(x, 1), 1)-0.5), corrected_y, gids);
for i = 1:4
    plot([i-linelength/2 i+linelength/2],...
        median(corrected_y(x ==i))*ones(1, 2), 'k');
end
xlabel('cell types');
ylabel('Shift time (s)');
xticks(1:5)
xticklabels({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});
ylim([-0.04 0.08]);
%% Evaluate the combination of both factors
figure;
subplot(1, 3, 1); hold on
x = Stretching_info(:, 1);
y = -tf_delay_time(:, 1);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('Stretching factor');
ylabel('Moving spot delay (s)');
title(sprintf('Shift & stretch correlation (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:)))); % Width(1, 3) Height (2)

subplot(1, 3, 2); hold on
x = - Stretching_info(delay_info(:, 8), 1).*tf_delay_time(delay_info(:, 8), 1);
y = delay_info(:, 1).*delay_info(:, 3);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('Stretching factor * Time shift');
ylabel('Moving spot delay distance (um)');
title(sprintf('Height difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:)))); % Width(1, 3) Height (2)

subplot(1, 3, 3); hold on
x = - Stretching_info(delay_info(:, 8), 1).*tf_delay_time(delay_info(:, 8), 1);
y = delay_info(:, 1);
[yFit, beta] = SimpleRegression2D(x, y);
[~, sids] = sort(x);
scatter(x, y, 25, 'k', 'filled');
plot(x(sids), yFit(sids), '--r');
xlabel('Stretching factor * Time shift');
ylabel('Moving spot delay (s)');
title(sprintf('Height difference (slope: %.03G, r=%.03G)', beta(2), corr(x(:), y(:)))); % Width(1, 3) Height (2)

sgtitle(sprintf('Direction %d', direc_ids(di)));


%% Visualize (corrected shift) responses according to temporal filter shift
cids = DataInfo(:, 3); % retina specific
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
colors = lines(length(ucids));
up_sample_fac = 10;
IsDisplay =0;
pixel_size = 2.5;
diame_exam_ids = 6:7;
speed_exam_ids = 3:5;
anchor_pos = [400 300];
max_cell_per_type = 5;
cell_type_list = [0, 1, 10, 11];
is_corrected = 1;
% close all
ii = 1;
switch direc_ids(di)
    case 0
        tfdelay_to_shiftlag = 1;
        shift_pos_ind = 1;
    case 90
        tfdelay_to_shiftlag = -1;
        shift_pos_ind = 2;
    case 180
        tfdelay_to_shiftlag = -1;
        shift_pos_ind = 1;
end
gather_responses = nan(max_cell_per_type*length(diame_exam_ids)*length(speed_exam_ids)*4, size(summary_trace, 2)*up_sample_fac);
gather_info = nan(max_cell_per_type*length(diame_exam_ids)*length(speed_exam_ids)*4, 3);
gi = 1;
for u = ucids(:)'
    clc
    fprintf('progress... %d/%d \n', u, ucids(end));
    for k =  1:length(diame_exam_ids) % length(diame_ids)-DisplaySquare:length(diame_ids)
        for j = 1:length(speed_exam_ids) % length(speed_ids)-DisplaySquare:length(speed_ids)
            cids = find(summary_info(:, 1) == di & summary_info(:, 4) == diame_exam_ids(k) &...
                summary_info(:, 5) == speed_exam_ids(j) & summary_info(:, 6) == u);
            type_ids = summary_info(cids, 3);
            ncell = length(cids);
            for p = 1:ncell
                % find the shift based on its temporal kernel shift
                % from the temporal kernel template of that type
                cell_id = summary_info(cids(p), 7); %
                s2 = summary_trace(cids(p), :);
                exist_time = all(~isnan(s2), 1);
                s2 = s2(exist_time);
                t = (0:length(s2)-1)/Fz;
                s2_up = SimpleSignalUpSample(s2, up_sample_fac);
                t_up = (0:length(s2_up)-1)/(Fz*up_sample_fac);

                offset_time = -tf_delay_time(cell_id, 1);
                iLag = round(tfdelay_to_shiftlag*Fz*up_sample_fac*offset_time + length(s2_up));
                s2_up_shift = circshift(s2_up, [0 iLag]);

                if is_corrected
                    gather_responses(gi, 1:length(s2_up)) = s2_up_shift;
                else
                    gather_responses(gi, 1:length(s2_up)) = s2_up;
                end
                gather_info(gi, :) = [type_ids(p), diame_exam_ids(k), speed_exam_ids(j)];
                gi = gi +1;
                %plot(t_up, s2_up, 'Color', line_color);
                % plot(t_up, s2_up_shift, 'Color', line_color);
                if IsDisplay
                    close all
                    figure; hold on
                    plot(t_up, s2_up, 'k');
                    plot(t_up, s2_up_shift, 'm');
                    title(sprintf('displacement: %3G (um)', rf_offset));
                    keyboard;
                end
            end
        end
    end
end


figure;
ii = 1;
for k =  1:length(diame_exam_ids) % length(diame_ids)-DisplaySquare:length(diame_ids)
    for j = 1:length(speed_exam_ids) % length(speed_ids)-DisplaySquare:length(speed_ids)
        subplot(length(diame_exam_ids), length(speed_exam_ids), ii); hold on
        ii = ii + 1;
        for q = 1:length(cell_type_list)
            switch cell_type_list(q)
                case 10 % ON-Nasal
                    line_color = [92 204 206]/255;
                case 11 % ON-Temporal
                    line_color = [255 139 139]/255;
                case 0 % OFF-Nasal
                    line_color = [60 91 89]/255;
                case 1 % OFF-Temporal
                    line_color = [136 74 75]/255;
            end
            cids = gather_info(:, 1) == cell_type_list(q) & gather_info(:, 2) == diame_exam_ids(k) &...
                gather_info(:, 3) == speed_exam_ids(j);
            avg = mean(gather_responses(cids, :), 1);
            sem = std(gather_responses(cids, :), [], 1)/sqrt(sum(cids));
            exist_time = all(~isnan([avg; sem]), 1);
            avg = avg(exist_time);
            sem = sem(exist_time);
            t = (0:length(avg)-1)/(Fz*up_sample_fac);
            shadePlot(t, avg, sem*1.96, line_color);
        end

        if k == 1
            title(sprintf('Speeds:%d(um/s)', speed_ids(speed_exam_ids(j))))
            Label_spds{j} = sprintf('%d (um/s)', speed_ids(speed_exam_ids(j)));
        end
        if j == 1
            if k ==1
                ylabel(sprintf('Diameters:%d (um)', diame_ids(diame_exam_ids(k))));
            else
                ylabel(sprintf('%d (um)', diame_ids(diame_exam_ids(k))));
            end
            Label_dias{k} = sprintf('%d (um)', diame_ids(diame_exam_ids(k)));
        elseif j == 2
            ylabel('Firing rate (spike/s)')
        end
    end
end
if is_corrected
    sgtitle(sprintf('Direction:%d (corrected)', direc_ids(di)));
else
    sgtitle(sprintf('Direction:%d (uncorrected)', direc_ids(di)));
end

%% [Section D] Calculate delay based on ON or OFF template
cids = DataInfo(:, 3); % retina specific
ucids = unique(cids);
ucids = ucids(~isnan(ucids));
colors = lines(length(ucids));
DisplaySquare = 1; % the last DisplaySquare + 1
up_sample_fac = 10;
IsDisplay = 0;
pixel_size = 2.5;
diame_exam_ids = 6:7;
speed_exam_ids = 3:5;
close all
di = 1 % 1:length(direc_ids)
ii = 1;
posi_info = [];
delay_info = [];

for k =  1:length(diame_exam_ids) % length(diame_ids)-DisplaySquare:length(diame_ids)
    for j = 1:length(speed_exam_ids) % length(speed_ids)-DisplaySquare:length(speed_ids)
        for q = 1:length(cell_type_list)
            cids = find(summary_info(:, 1) == di & summary_info(:, 4) == diame_exam_ids(k) &...
                summary_info(:, 5) == speed_exam_ids(j) & ismember(summary_info(:, 3), cell_type_list(q)));
            ncell = length(cids);

            % switch cell_type_list(q)
            %     case {0, 1}
            %         contrast_type = [1, 2];
            %     case {10, 11}
            %         contrast_type = [3, 4];
            % end

            tids = template_info(:, 1) == di & ismember(template_info(:, 2), q) &...
                template_info(:, 3) == diame_exam_ids(k) & template_info(:, 4) == speed_exam_ids(j);
            for p = 1:ncell
                s1 = mean(template_traces(tids,  :), 1);
                s2 = summary_trace(cids(p), :);
                exist_time = all(~isnan([s1; s2]), 1);
                s1 = s1(exist_time);
                s2 = s2(exist_time);
                t = (0:length(s1)-1)/Fz;

                s1_up = SimpleSignalUpSample(s1(exist_time), up_sample_fac);
                s2_up = SimpleSignalUpSample(s2(exist_time), up_sample_fac);
                t_up = (0:length(s1_up)-1)/(Fz*up_sample_fac);

                [c, lags] = xcorr(s1_up(:), s2_up(:));
                [~, iLag] = max(c);
                s3_up = circshift(s2_up, [0 iLag]);
                shift_time = lags(iLag)/(Fz*up_sample_fac);

                pos_1 = [400 300];
                pos_2 = cent_cells(summary_info(cids(p), 7), :);
                posi_info = [posi_info; pos_2, pos_1];
                delay_info = [delay_info; shift_time, diame_ids(diame_exam_ids(k)), speed_ids(speed_exam_ids(j)),...
                    diame_exam_ids(k), speed_exam_ids(j), q, cell_type_list(q), summary_info(cids(p), 7)];

                if IsDisplay
                    close all
                    figure; hold on
                    plot(t, s1, '--k');
                    plot(t, s2, '--m');
                    plot(t_up, s1_up, 'k');
                    plot(t_up, s2_up, 'm');
                    plot(t_up, s3_up, 'b');
                    title(sprintf('Lag: %.03G (s)', shift_time));
                    keyboard;
                end
            end

        end
    end
end


%% Calculate TF kernel delay time to template of ON and OFF cells
exam_measurement_id = 'stretching'; 
switch lower(exam_measurement_id)
    case 'shifting'
        column_id = 6;
        y_label = 'TF delay' ;
    case 'stretching'
        column_id = 4;
        y_label = 'TF stretch' ;
end

close all
x = Stretching_info(delay_info(:, 8), column_id);
y = -delay_info(:, 1).*delay_info(:, 3);
z = (cell_type_list == delay_info(:, 7))*(1:4)';

Colors = nan(4, 3);
for i = 1:4
    switch i
        case 3 % ON-Nasal % 10
            line_color = [92 204 206]/255;
        case 4 % ON-Temporal % 11
            line_color = [255 139 139]/255;
        case 1 % OFF-Nasal %  0
            line_color = [60 91 89]/255;
        case 2 % OFF-Temporal % 1
            line_color = [136 74 75]/255; % [136 74 75]
    end
    Colors(i, :) = line_color;
end

figure; 
subplot(1, 3, 1);
gscatter(x, y, z, Colors, [], 15);
box off
xlabel(sprintf('Avg. %s', y_label));
ylabel('Avg. moving spot delay');
legend({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});

unique_cell_ids = unique(delay_info(:, 8));
num_unique_cell_ids = numel(unique_cell_ids);
x = nan(num_unique_cell_ids, 1);
y = nan(num_unique_cell_ids, 1);
z = nan(num_unique_cell_ids, 1);

for i = 1:num_unique_cell_ids
    uids = delay_info(:, 8) == unique_cell_ids(i);
    x(i) = mean(Stretching_info(delay_info(uids, 8), column_id)); % 6
    y(i) = mean(-delay_info(uids, 1).*delay_info(uids, 3));
    z(i) = find(cell_type_list == mean(delay_info(uids, 7)));

end


subplot(1, 3, 2); hold on
gscatter(x, y, z, Colors, [], 25);
corr_off_ids = ismember(z, 1:2);
[~, beta] = SimpleRegression2D(x(corr_off_ids), y(corr_off_ids));
[yFit, ~] = SimpleRegression2D(x, [], beta);
[~, sids] = sort(x);
plot(x(sids), yFit(sids), '--k');
box off
xlabel(sprintf('Avg. %s', y_label));
ylabel('Avg. moving spot delay');
legend({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});

subplot(1, 3, 3); hold on
gscatter(x, y-yFit, z, Colors, [], 25);
plot(x(sids), yFit(sids)-yFit(sids), '--k');
box off
legend({'OFFnasal', 'OFFtemporal', 'ONnasal', 'ONtemporal'});
xlabel(sprintf('Avg. %s', y_label));
ylabel('Avg. moving spot delay');

