clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\StationarySpot\';
recording_id = 47;

%% ON - Nasal 51?, 49, Temporal (47, 46, 44), [43 42]
%% OFF - 50
switch recording_id
    case 1
        recordingname = 'a120723';
        load([ephys_data_folder recordingname '_0010_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_StationalSpot.mat']);
    case 2 % not yet finished
        recordingname = 'b120723';
        load([ephys_data_folder recordingname '_0002_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_StationalSpot.mat']);
        MinPeakHeight = 10;
        Is_extracellular = 0;
    case 4
        recordingname = 'd010224';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_2_StationalSpot.mat']);
    case 5 % far from center
        recordingname = 'a030124';
        load([ephys_data_folder recordingname '_0003_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 6 % far from center
        recordingname = 'c030124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 400;
    case 10
        recordingname = 'a030624';
        load([ephys_data_folder recordingname '_0002_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 11
        recordingname = 'a030924';
        load([ephys_data_folder recordingname '_0000_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 100;
    case 12
        recordingname = 'b030924';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 100;
    case 13
        recordingname = 'c030924';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 100;
    case 43
        recordingname = 'd030924';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 250; 
    case 14
        recordingname = 'a032924';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 150;
    case 21
        recordingname = 'b032924';
        load([ephys_data_folder recordingname '_0004_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_004_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 23
        recordingname = 'c032924';
        load([ephys_data_folder recordingname '_0003_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 400;
    case 22
        recordingname = 'c032924';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 1000;
    case 15
        recordingname = 'a040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 300;
    case 16
        recordingname = 'b040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 1500;
    case 25
        recordingname = 'c040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 1000;
    case 26
        recordingname = 'd040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 1500;
    case 27
        recordingname = 'e040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 2000;
    case 28
        recordingname = 'f040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 2000;
    case 29
        recordingname = 'g040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 2000;
    case 17
        recordingname = 'a040224';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 18
        recordingname = 'b040224';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 19
        recordingname = 'c040224';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 150;
    case 20
        recordingname = 'd040224';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 24
        recordingname = 'e040224';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 30
        recordingname = 'a040524';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 31
        recordingname = 'b040524';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 32
        recordingname = 'c040524';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 33
        recordingname = 'd040524';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 210;
    case 34
        recordingname = 'e040524';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 35
        recordingname = 'a040624';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 100;
    case 36
        recordingname = 'b040624';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 37 % probably dead
        recordingname = 'c040624';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 200;
    case 38
        recordingname = 'c040624';
        load([ephys_data_folder recordingname '_0002_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 400;
    case 39
        recordingname = 'd040624';
        load([ephys_data_folder recordingname '_0002_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 700;
    case 40
        recordingname = 'e040624';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 41
        recordingname = 'f040624';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 300;
    case 42
        recordingname = 'g040624';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 44
        recordingname = 'a081024';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 250;
    case 45
        recordingname = 'b081024';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 250;
    case 46
        recordingname = 'c081024';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 300;
    case 47
        recordingname = 'd081024';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 48
        recordingname = 'a081224';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 49
        recordingname = 'c081224';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 50
        recordingname = 'e081224';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
    case 51
        recordingname = 'f081224';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
end

%%
OUT = SS_OUT;
IN = SS_IN;
%OLED.height = 600;
%OLED.width = 800;
%OLED.pixelSize = 3.5; %size of each pixel in micron on retina
clear SS_OUT SS_IN
assert(abs(length(sectionData)/(OUT.EndTime-OUT.StartTime) -10000)<1);

%% Detrend
close all
figure; 
subplot(2, 1, 1); hold on
plot(sectionData, 'k');
% plot(smoothdata(sectionData, 'movmedian', 1*10000), 'r');
lpass = lowpass(sectionData, 0.05, 10000, 'Steepness', 0.99);
plot(lpass, 'r');

subplot(2, 1, 2); hold on
plot(sectionData, 'k');
plot(sectionData-lpass, 'r');
%%
% keyboard;
%%
sectionData = sectionData-lpass;

%%
Fz = 100;
WinT = [-0.5 0];
% convert trace to spike or firing rate
figure; hold on
if ~exist('Is_extracellular', 'var')
    Is_extracellular = 1;
end
if Is_extracellular
    [pks, locs] = findpeaks(-sectionData, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 40);
    plot(sectionData(1:locs(end)), 'k');
    plot(locs(1:end), -pks(1:end), 'rx');

else
    [pks, locs] = findpeaks(sectionData, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 40);
    plot(sectionData(1:locs(1000)), 'k');
    plot(locs(1:1000), pks(1:1000), 'rx');
end

%%
fs = 10000;
signals = zeros(size(sectionData));
signals(locs) = 1;
% t = (0:length(sectionData)-1)/fs;
downsample_size = round(fs/Fz);
add_more = mod(length(signals), downsample_size);
if add_more ~= 0
    add_more = downsample_size - add_more;
    sig = [signals(:)' zeros(1, add_more)];
else
    sig = signals(:)';
end
assert(mod(length(sig), downsample_size) == 0);
sig = reshape(sig(:), downsample_size, []);
sig = sum(sig, 1);
t = (0:length(sig)-1)/Fz;
ssigs = smoothdata(sig, 'gaussian', 0.05, 'SamplePoints', t');
figure; hold on
plot(t, sig, 'k');
plot(t, ssigs, 'b');
sig = ssigs*Fz;
clear ssigs

%%
%keyboard;
%% get new table with time ids
% [      1      ,      2     ,       3     ,      4     ,     5   ,    6,  7 ]
% [Global repeat, Stimulus id, Local repeat, contrast id, Diameter, Contrast, Time ]
OUT.DataTable(:, 8) = round(Fz*(OUT.DataTable(:, 7)-OUT.StartTime))+1;
assert(((OUT.EndTime-OUT.StartTime)*Fz-max(OUT.DataTable(:, 8)))>=IN.ISI*Fz);
diame_ids = unique(OUT.DataTable(:, 5));
contr_ids = unique(OUT.DataTable(:, 4));
Dtab = nan(length(diame_ids)*length(contr_ids), 4);
ii = 1;
tids = [2 1];
sig_add_time = round(IN.presentT(tids)*Fz);
num_repeat = IN.globalrepeat*IN.localrepeat;
for i = 1:length(diame_ids)
    for j = 1:length(contr_ids)
        ctime = OUT.DataTable(OUT.DataTable(:, 5) == diame_ids(i) & OUT.DataTable(:, 4) == contr_ids(j), 8);
        assert(length(ctime) == num_repeat);
        for v = 1:num_repeat
            Dtab(ii, :) = [i, j, ctime(v), ctime(v)+sig_add_time(j)];
            ii = ii + 1;
        end
    end
end
Dtab(:, 5) = (Dtab(:, 4)-Dtab(:, 3))/Fz;



%% only global repeat trace
MaxTwin = sum(IN.presentT)*IN.localrepeat+1; % in second (3.2 the max stimulus presentation length)
TimeBefore = 0.5;
nt = round((MaxTwin+TimeBefore)*Fz)+1;
ct = (0:nt-1)/Fz;
DataG = nan(length(diame_ids), IN.globalrepeat, nt);
figure;
for i = 1:length(diame_ids)
    subplot(1, length(diame_ids),  i); hold on
    cids = find(Dtab(:, 1)==i & Dtab(:, 2)==2);
    gsigs = nan(IN.globalrepeat, nt);
    for v = 1:IN.globalrepeat
        ctim = Dtab(cids((v-1)*IN.localrepeat+1), 3);
        sids = ctim(1)-TimeBefore*Fz:ctim(1)+round(MaxTwin*Fz);
        if sids(end) > length(sig)
            sids = sids(1):length(sig);
        end
        csig = sig(sids);
        gsigs(v, 1:length(sids)) = csig';

        plot(ct(1:length(sids)), csig, 'Color', 0.5*ones(1, 3));
        xlim([0, ct(end)]);
    end
    DataG(i, :, :) = gsigs;
end
%%
MaxTwin = max(IN.presentT)+0.3; % in second (3.2 the max stimulus presentation length)
TimeBefore = 0.5;
%%

% close all
figure;

nt = round((MaxTwin+TimeBefore)*Fz)+1;
ct = (0:nt-1)/Fz;
contrast_name = {'OFF', 'ON'};
Data = nan(length(diame_ids), length(contr_ids), num_repeat-1, nt);
for j = 1:length(contr_ids)
    for i = 1:length(diame_ids)
        subplot(length(contr_ids), length(diame_ids),  (j-1)*length(diame_ids)+i); hold on
        cids = find(Dtab(:, 1)==i & Dtab(:, 2)==j);
        % if j == 2
        %     cids(1) = [];
        % elseif j == 1
        %     cids(end) = [];
        % end
         
        gsigs = nan(num_repeat-1, nt);
        for v = 1:num_repeat-1
            ctim = Dtab(cids(v), 3);
            csig = sig(ctim(1)-TimeBefore*Fz:ctim(1)+round(MaxTwin*Fz));
            gsigs(v, :) = csig';
            
            plot(ct, csig, 'Color', 0.5*ones(1, 3));
            xlim([0, ct(end)]);
        end
        Data(i, j, :, :) = gsigs;
        plot(ct, mean(gsigs, 1), 'k');
        plot([0 0.5], 200*ones(1, 2), 'k', 'LineWidth', 3);
        plot([0.5 2], 200*ones(1, 2), 'y', 'LineWidth', 3);
        plot([2 ct(end)], 200*ones(1, 2), 'k', 'LineWidth', 3);
        if i == 1
            ylabel(sprintf('%s', contrast_name{j}));
        end
        if j == length(contr_ids)
            xlabel('Time (s)');
        end
        if j == 1
            if i == 1
                title(sprintf('Diameters:%d (um)', IN.diameters(i)));
            else
                title(sprintf('%d (um)', IN.diameters(i)));
            end
        end
        if i == 2
            ylabel('Firing rate (spike/s)')
        end
        ylim([0 200]);
    end
end
sgtitle(sprintf('%s Stationary Spot', recordingname));
%%
diameters = IN.diameters;
ISI = IN.ISI;
saveFileName = sprintf('%s_stationary_spot.mat', recordingname);
save(sprintf('%s/%s', save_data_folder, saveFileName), 'diameters', 'Fz', 'nt','ct', 'contrast_name',...
    'num_repeat', 'Dtab', 'ISI', 'Data', 'DataG', 'contr_ids', 'diame_ids', 'num_repeat', 'recordingname');
%% Analysis ends here
keyboard;
%%
Win_LON = [0.2 0.5];
Win_LOFF = [0.7 1];
cData = nan(length(diame_ids), length(contr_ids), IN.repeat, 2);
Data = nan(length(diame_ids), length(contr_ids), IN.repeat, round(MaxTwin*Fz));
Moving_end_time_ids = nan(length(diame_ids), length(contr_ids), IN.repeat);
close all

figure;
for i = 1:length(diame_ids)
    for j = 1:length(contr_ids)
        subplot(length(diame_ids), length(contr_ids), (i-1)*length(contr_ids)+j); hold on
        cids = find(Dtab(:, 1)==i & Dtab(:, 2)==j);
        gsigs = nan(IN.repeat, round(MaxTwin*Fz));
        for v = 1:IN.repeat
            ctim = Dtab(cids(v), 3:4);
            csig = sig(ctim(1):ctim(1)+round(MaxTwin*Fz));
            if length(csig) > round(MaxTwin*Fz)
                gsigs(v, :) = csig(1:round(MaxTwin*Fz))';
            else
                gsigs(v, 1:length(csig)) = csig(:)';
            end
            ct = (0:length(csig)-1)/Fz;
            plot(ct, csig, 'Color', 0.5*ones(1, 3));
            cEndTid = ctim(2)-ctim(1)+1;
            Moving_end_time_ids(i, j, v) = cEndTid;
            plot(ct(cEndTid)*ones(1, 2), [0 60], 'r');
            xlim([0, MaxTwin]);
            
            rid = ct>Win_LON(1) & ct<=Win_LON(2);
            cData(i, j, v, 1) = mean(csig(rid));
            rid = ct>Win_LOFF(1) & ct<=Win_LOFF(2);
            cData(i, j, v, 2) = mean(csig(rid));
        end
        Data(i, j, :, :) = gsigs;
        plot(ct(1:round(MaxTwin*Fz)), mean(gsigs, 1), 'k');
        if i == 1
            title(sprintf('Contrasts:%2G(%%)', round(100*(IN.contrasts(j)-0.5))))
        elseif i == length(diame_ids)
            xlabel('Time (s)');
        end
        if j == 1
            if i == 1
                ylabel(sprintf('Diameters:%d (um)', IN.diameters(i)));
            else
                ylabel(sprintf('%d (um)', IN.diameters(i)));
            end
        elseif j == 2
            ylabel('Firing rate (spike/s)')
        end
        ylim([0 200]);
    end
end
%%
dim1_diameters = IN.diameters;
dim2_contrasts = IN.contrasts;
dim3_repeats = 1:IN.repeat;
dim4_time = (0:size(Data, 4)-1)/Fz;
saveFileName = sprintf('%s_moving_stationary_spot_contrast_processed.mat', recordingname);
save(sprintf('./ProcessedData/%s', saveFileName), 'dim1_diameters', 'dim2_contrasts', 'dim3_repeats',...
'dim4_time', 'Data', 'Moving_end_time_ids');
%% Contrast steps
cts = round(100*(IN.contrasts-0.5));

figure;
for i = 1:length(diame_ids)
    subplot(2, 4, i); hold on
    mc = squeeze(mean(cData(i, :, :, 1), 3));
    stdc = squeeze(std(cData(i, :, :, 1), [], 3))/sqrt(IN.repeat);
    errorbar(cts, mc, stdc, 'Color', 0.5*ones(1, 3), 'CapSize', 0);
    mc = squeeze(mean(cData(i, :, :, 2), 3));
    stdc = squeeze(std(cData(i, :, :, 2), [], 3))/sqrt(IN.repeat);
    errorbar(cts, mc, stdc, 'k', 'CapSize', 0);
    box off
    ylim([0 150]);
    xlabel('Contrast (%)');
    ylabel('Average firing rate (spike/s)');
    if i == 1
        legend({'Light ON', 'Light OFF'});
    end
    title(sprintf('Spot diameter:%d (um)',IN.diameters(i)));
end
