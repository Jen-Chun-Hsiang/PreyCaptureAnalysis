clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_name = 'a081024';
load_recording_name = [recording_name '01'];
Fz = 100;
%% get the linear and nonlinear parts
load([load_data_folder load_recording_name '.mat'], 'masked_STAmat', 'PBs', 'FRs');
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_MovingNoise_1.mat'], 'OLED');
% Spatial temporal filter is in "masked_STAmat"

% Nonlinear parts can be derived from "PBs", "FRs"

[nl_fuc, divider]  = getNonlinearFunc(PBs, FRs);
%% Simulation (0) Step full field

x = -1*ones(1, Fz*2);
x(51:150) = 1;
cmov = -1*ones(size(masked_STAmat));
num_time = length(x);
t = (0:num_time-1)/Fz;
resp = nan(1, num_time);
dim1 = size(cmov, 1);
dim2 = size(cmov, 2);
for i = 1:num_time
    resp(i) = mean(cmov.*masked_STAmat, 'all');
    cmov(:, :, 1) = [];
    cmov(:, :, end+1) = x(i)*ones(dim1, dim2);
    clc
    fprintf('Progress ... %d/%d \n', i, num_time);
end
%%
figure;
subplot(1, 2, 1)
plot(t, resp, 'k'); box off
ylabel('linear effective contrast (arbi. unit');
subplot(1, 2, 2)
plot(t, nl_fuc(resp/divider), 'k'); box off
ylabel('Predicted firing rate (spike/s)')

%% Simulation (1) varied size spots
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_VariedSizeSpot.mat'], 'SS_IN');
spot_size = SS_IN.diameters;
pix2um = OLED.pixelSize;
[X, Y] = meshgrid(1:size(masked_STAmat, 2), 1:size(masked_STAmat, 1));
masked = sqrt((X-0.5*(size(masked_STAmat, 2)+1)).^2 +  (Y-0.5*(size(masked_STAmat, 1)+1)).^2) < spot_size(6)/pix2um;
%%
figure;
subplot(1, 2, 1)
imagesc(masked);
subplot(1, 2, 2)
imagesc(masked_STAmat(:, :, 45));

%% load processed data of responses
loadFileName = sprintf('%s_stationary_spot.mat', recording_name);
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\StationarySpot\';
load(sprintf('%s/%s', load_data_folder, loadFileName), 'diameters','ct', 'contrast_name',...
    'num_repeat', 'ISI', 'Data', 'DataG', 'contr_ids', 'diame_ids');
% Data [diameters, contrasts, repeat, time]
%%
dia_id = 5;
figure;
plot(ct, squeeze(mean(Data(dia_id, 2, :, :), 3, 'omitnan')));
title(sprintf('%d (um)', diameters(dia_id)));

%%
figure; hold on
x = mean(Data(:, 2, :, ct>0.1 & ct<0.5), 4, 'omitnan');
plot(diameters, squeeze(mean(x, 3, 'omitnan')), 'Color', [255 210 0]/255);
x = mean(Data(:, 1, :, ct>0.1 & ct<0.5), 4, 'omitnan');
plot(diameters, squeeze(mean(x, 3, 'omitnan')), 'k');
ylabel('Firing rate (spike/s)');
xlabel('Time (s)');
legend({'ON', 'OFF'});

%%
figure; hold on
x = mean(Data(:, 2, :, ct>0.1 & ct<0.5), 4, 'omitnan') - mean(Data(:, 2, :, ct<0.05), 4, 'omitnan');
plot(diameters, squeeze(mean(x, 3, 'omitnan')), 'Color', [255 210 0]/255);
x = mean(Data(:, 1, :, ct>0.1 & ct<0.5), 4, 'omitnan') - mean(Data(:, 1, :, ct<0.05), 4, 'omitnan');
plot(diameters, squeeze(mean(x, 3, 'omitnan')), 'k');
ylabel('delta Firing rate (spike/s)');
xlabel('Time (s)');
legend({'ON', 'OFF'});

%%
keyboard;
%%
gray_mask = sqrt((X-0.5*(size(masked_STAmat, 2)+1)).^2 +  (Y-0.5*(size(masked_STAmat, 1)+1)).^2) < 0.5*SS_IN.MaskSize/pix2um;
close all
for j = 1:length(diame_ids)
    masked = sqrt((X-0.5*(size(masked_STAmat, 2)+1)).^2 +  (Y-0.5*(size(masked_STAmat, 1)+1)).^2) < 0.5*diameters(j)/pix2um;
    x = zeros(1, Fz*13.5);
    x(51:200) = 1;
    x(201:350) = -1;
    x(351:500) = 1;
    x(501:650) = -1;
    x(651:800) = 1;
    x(801:950) = -1;
    x(951:1100) = 1;
    x(1101:1250) = -1;
    canvas = -1*ones(size(masked));
    canvas(gray_mask == 1) = 0;
    cmov = repmat(canvas, 1, 1, size(masked_STAmat, 3));
    num_time = length(x);
    t = (0:num_time-1)/Fz;
    if j == 1
        resp_ON1 = nan(length(diame_ids), num_time);
    end
    dim1 = size(cmov, 1);
    dim2 = size(cmov, 2);
    for i = 1:num_time
        resp_ON1(j, i) = mean(cmov.*masked_STAmat, 'all');
        cmov(:, :, 1) = [];
        canvas = -1*ones(size(masked));
        canvas(gray_mask == 1) = 0;
        canvas(masked == 1) = x(i);
        cmov(:, :, end+1) = canvas;
        % if mod(i, 10) == 1
        %     figure(1);
        %     imagesc(canvas, [-1 1]); colorbar;
        % end
        clc
        fprintf('Progress ON1... %d/%d, %d/%d \n', j, length(diame_ids), i, num_time);
    end
end
for j = 1:length(diame_ids)
    masked = sqrt((X-0.5*(size(masked_STAmat, 2)+1)).^2 +  (Y-0.5*(size(masked_STAmat, 1)+1)).^2) < 0.5*diameters(j)/pix2um;
    x = zeros(1, Fz*13.5);
    x(51:200) = -1;
    x(201:350) = 1;
    x(351:500) = -1;
    x(501:650) = 1;
    x(651:800) = -1;
    x(801:950) = 1;
    x(951:1100) = -1;
    x(1101:1250) = 1;
    canvas = -1*ones(size(masked));
    canvas(gray_mask == 1) = 0;
    cmov = repmat(canvas, 1, 1, size(masked_STAmat, 3));
    num_time = length(x);
    t = (0:num_time-1)/Fz;
    if j == 1
        resp_OFF1 = nan(length(diame_ids), num_time);
    end
    dim1 = size(cmov, 1);
    dim2 = size(cmov, 2);
    for i = 1:num_time
        resp_OFF1(j, i) = mean(cmov.*masked_STAmat, 'all');
        cmov(:, :, 1) = [];
        canvas = -1*ones(size(masked));
        canvas(gray_mask == 1) = 0;
        canvas(masked == 1) = x(i);
        cmov(:, :, end+1) = canvas;
        % if mod(i, 10) == 1
        %     figure(1);
        %     imagesc(canvas, [-1 1]); colorbar;
        % end
        clc
        fprintf('Progress OFF1... %d/%d, %d/%d \n', j, length(diame_ids), i, num_time);
    end
end
%%
save_file_name = sprintf('%s_%s_stationaryspot_simulated.mat', recording_name);

%%
% save(sprintf('./Results/MovingBar/%s', save_file_name), 'diame_ids', 'x',...
%     'resp_ON1', 'resp_OFF1');
%%
keyboard;
%% 
load(sprintf('./Results/MovingBar/%s', save_file_name), 'diame_ids', 'x',...
    'resp_ON1', 'resp_OFF1');
%%
colors = parula (4);
figure; 
for i = 1:4
    subplot(1, 2, 1); hold on
    plot(squeeze(Data(7, 2, i, :)), 'Color', colors(i, :));
    subplot(1, 2, 2); hold on
    plot(squeeze(Data(7, 2, i+4, :)), 'Color', colors(i, :));
end
%%
c_dia_id = 5;
figure; 
subplot(1, 2, 1);hold on
a = squeeze(DataG(c_dia_id, 1, 151:end));
plot(a, 'Color', 'k');
b = squeeze(resp_ON1(c_dia_id, 151:end));
plot(b*max(a)./max(b), 'Color','b');

subplot(1, 2, 2);hold on
a = squeeze(DataG(c_dia_id, 2, 1:end-150));
plot(a, 'Color','k');
b = squeeze(resp_OFF1(c_dia_id, 151:end));
plot(b*max(a)./max(b), 'Color', 'b');

%%
figure; 
%%
figure;
plot(nl_fuc(resp(1, :)/divider));
%% process simulation data
nt = size(Data, 4);
Data_s = nan(length(diame_ids), 2, nt);
for i = 1:length(diame_ids)
    Data_s(i, 2, :) = resp(i, 301:300+nt);
    Data_s(i, 1, :) = resp(i, 151:150+nt);
end

%%
cData_s = Data_s*150./max(Data_s(:));
ct = (0:(nt-1))/Fz;
figure;
for i = 1:2
    for j = 1:length(diame_ids)
        subplot(2, length(diame_ids), (i-1)*length(diame_ids) + j);
        hold on
        if i == 2
            plot(ct, squeeze(mean(Data(j, 2, :, :), 3, 'omitnan')), 'k');
            plot(ct, nl_fuc(squeeze(Data_s(j, 2, :))/divider), 'r');
            plot(ct, 0.6*squeeze(cData_s(j, 2, :))+50, 'b');
            xlabel('Time (s)');
        elseif i == 1
            plot(ct, squeeze(mean(Data(j, 1, :, :), 3, 'omitnan')), 'k');
            plot(ct, nl_fuc(squeeze(Data_s(j, 1, :))/divider), 'r');
            plot(ct, 0.6*squeeze(cData_s(j,1, :))+50, 'b');
            title(sprintf('%d (um)', diameters(j)));
        end
        ylim([0 150]);
        if j == 1
            if i == 1
                ylabel(sprintf('OFF \n Firing rate (spike/s)'));
            elseif i == 2
                ylabel(sprintf('ON \n Firing rate (spike/s)'));
            end
        end
    end
end
sgtitle(sprintf('%s', recording_name));

%% 