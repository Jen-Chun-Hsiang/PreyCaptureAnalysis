clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\StationarySpot\';
recording_id = 15;
switch recording_id
    case 1
        recordingname = 'a120723';
        load([ephys_data_folder recordingname '_0010_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_StationalSpot.mat']);
    case 4
        recordingname = 'd010224';
        load([ephys_data_folder recordingname '_0001_3.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_2_StationalSpot.mat']);
    case 15
        recordingname = 'a040124';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_VariedSizeSpot.mat']);
        MinPeakHeight = 500;
end

%%
OUT = SS_OUT;
IN = SS_IN;
OLED.height = 600;
OLED.width = 800;
OLED.pixelSize = 3.5; %size of each pixel in micron on retina
clear SS_OUT SS_IN
assert(abs(length(sectionData)/(OUT.EndTime-OUT.StartTime) -10000)<1);

%%
Fz = 100;
WinT = [-0.5 0];
% convert trace to spike or firing rate
[pks, locs] = findpeaks(-sectionData, 'MinPeakHeight', 0, 'MinPeakDistance', 40);
figure; hold on
% plot(sectionData(1:locs(1000)), 'k');
% plot(locs(1:1000), pks(1:1000), 'rx');
plot(sectionData(1:locs(end)), 'k');
plot(locs(1:end), -pks(1:end), 'rx');
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
keyboard;
%% get new table with time ids
% [Direction, Speed, Diameter]
OUT.DataTable(:, 5) = round(Fz*(OUT.DataTable(:, 4)-OUT.StartTime))+1;
assert(((OUT.EndTime-OUT.StartTime)*Fz-max(OUT.DataTable(:, 5)))>=IN.ISI*Fz);
diame_ids = unique(OUT.DataTable(:, 2));
contr_ids = unique(OUT.DataTable(:, 3));
Dtab = nan(length(diame_ids)*length(contr_ids), 4);
ii = 1;
for i = 1:length(diame_ids)
    for j = 1:length(contr_ids)
        ctime = OUT.DataTable(OUT.DataTable(:, 2) == diame_ids(i) & OUT.DataTable(:, 3) == contr_ids(j), 5);
        assert(length(ctime) == 2*IN.repeat);
        for v = 1:IN.repeat
            Dtab(ii, :) = [i, j, ctime((v-1)*2+1), ctime(v*2)];
            ii = ii + 1;
        end
    end
end
Dtab(:, 5) = (Dtab(:, 4)-Dtab(:, 3))/Fz;

%%
MaxTwin = max(Dtab(:, 5))+IN.ISI-0.05; % in second (3.2 the max stimulus presentation length)
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
