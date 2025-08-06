clear; close all; clc;
DataFileFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily\Data\Stimulation';
DataFileName = 'Temporal_AlphaRGC_y030924_000_Retina_3_MovingSpot.mat';
Ds = load(fullfile(DataFileFolder, DataFileName));
DataFileFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\TestStimulus\MatFile';
DataFileName = 'TestHighSpeedMovingSpot_030824003.mat';
Dr = load(fullfile(DataFileFolder, DataFileName));
ImageData = Dr.ImageData;
clear Dr
%% Link trigger
channel_id = 3;
Sampling_rate = 5.92;
detect_threshold = 2e4;
minSectionDuration = 5;


Fz = size(ImageData, 1)*size(ImageData, 1)*Sampling_rate;
figure; 
Trigger = [];
for i = 1:size(ImageData, 3)
    Trigger = [Trigger; reshape(squeeze(ImageData(:, :, i, channel_id))', [], 1)];
end
t = (0:length(Trigger)-1)/Fz;
plot(t, Trigger);
minSamples = minSectionDuration * Fz; % Convert duration to samples
trigger_edges = tiggerdetection_continous(Trigger, detect_threshold);
trigger_edges(diff(trigger_edges, [], 2)<minSamples, :) = [];

channel_id = 1;
Signal = [];
for i = 1:size(ImageData, 3)
    Signal = [Signal; reshape(squeeze(ImageData(:, :, i, channel_id))', [], 1)];
end
figure; 
t = (0:length(Signal)-1)/Fz;
plot(t, Signal);
%%
DataTable = Ds.MS_OUT.DataTable;
DataTable = DataTable(ismember(DataTable(:, 1), [1, 2]), [1 3 4 5]);
DataTable(:, end+1) = round(trigger_edges(1) +...
    (trigger_edges(2)-trigger_edges(1))*(DataTable(:, end)-Ds.MS_OUT.StartTime)/(Ds.MS_OUT.EndTime - Ds.MS_OUT.StartTime));
%%
for i = 1:2 % different size
    cDataTable = DataTable(DataTable(:, 2) == i, :);
    nRepeat = size(cDataTable, 1)/2;
    Dur_length = cDataTable(2, 5) - cDataTable(1, 5);
    Resp = nan(nRepeat, Dur_length);
    cids = find(cDataTable(:, 1) == 1);
    cids_2 = find(cDataTable(:, 1) == 2);
    for j = 1:nRepeat
       sid = cDataTable(cids(j), 5);
       sid_2 = cDataTable(cids_2(j), 5);
       a = Signal(sid:sid_2-1);
       Resp(j, :) =  interp1(1:length(a), a,1:Dur_length);
    end
    figure; 
    plot(mean(Resp(2, :), 1));
    keyboard;
    
end

%%
[c, lags] = xcorr(Resp(3, :)', Resp(5, :)', 'unbiased');
figure; 
plot(lags, c);