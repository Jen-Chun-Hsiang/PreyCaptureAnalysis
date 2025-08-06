clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar\';

recording_id = 1;
switch recording_id
    case 1
        recordingname = 'a081024';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingBar_1.mat']);
        MinPeakHeight = 450;
    case 2
        recordingname = 'b081024';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingBar_1.mat']);
        MinPeakHeight = 300;
    case 3
        recordingname = 'c081024';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingBar_1.mat']);
        MinPeakHeight = 500;
    case 4
        recordingname = 'd081024';
        load([ephys_data_folder recordingname '_0001_2.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingBar_1.mat']);
        MinPeakHeight = 1000;

end

%%
OUT = MB_OUT;
IN = MB_IN;
OLED.height = 600;
OLED.width = 800;
OLED.pixelSize = 3.5; %size of each pixel in micron on retina
clear MS_OUT MS_IN
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
keyboard;

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
OUT.FrmTable(:, 7:9) = OUT.patternTable(OUT.FrmTable(:, 4), :);
OUT.FrmTable(:, 10) = round(Fz*(OUT.FrmTable(:, 5)-OUT.StartTime))+1;

%%
direc_ids = unique(OUT.FrmTable(:, 7));
speed_ids = unique(OUT.FrmTable(:, 8));
diame_ids = unique(OUT.FrmTable(:, 9));
Dtab = nan(length(direc_ids)*length(speed_ids)*length(diame_ids), 5);
ii = 1;
for i = 1:length(direc_ids)
    for j = 1:length(speed_ids)
        for k = 1:length(diame_ids)
            ctime = OUT.FrmTable(OUT.FrmTable(:, 7) == direc_ids(i) & OUT.FrmTable(:, 8) == speed_ids(j) &...
                OUT.FrmTable(:, 9) == diame_ids(k), 10);
            assert(length(ctime) == 2*IN.nRepeat);
            for v = 1:IN.nRepeat
                Dtab(ii, :) = [i, j, k, ctime((v-1)*2+1), ctime(v*2)];
                ii = ii + 1;
            end
        end
    end
end
Dtab(:, 6) = (Dtab(:, 5)-Dtab(:, 4))/Fz;
%%
MaxTwin = max(Dtab(:, 6))+IN.ISI; % in second (3.2 the max stimulus presentation length)

%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
cData = nan(length(direc_ids), length(diame_ids), length(speed_ids));
Data = nan(length(direc_ids), length(diame_ids), length(speed_ids), IN.nRepeat, round(MaxTwin*Fz));
Moving_end_time_ids = nan(length(direc_ids), length(diame_ids), length(speed_ids), IN.nRepeat);
clear Label_dias Label_spds
close all
for i = 1:length(direc_ids)
    figure;
    for k = 1:length(diame_ids)
        for j = 1:length(speed_ids)
            subplot(length(diame_ids), length(speed_ids), (k-1)*length(speed_ids)+j); hold on
            cids = find(Dtab(:, 1)==i & Dtab(:, 3)==k & Dtab(:, 2)==j);
            gsigs = nan(IN.nRepeat, round(MaxTwin*Fz));
            for v = 1:IN.nRepeat
                ctim = Dtab(cids(v), 4:5);
                tids = ctim(1):ctim(2)+round(IN.ISI*Fz);
                if max(tids) > length(sig)
                    csig = nan(1, length(tids));
                    csig(1:(length(sig)-ctim(1))+1) = sig(ctim(1):length(sig));
                else
                    csig = sig(tids);
                end
                if length(csig) > round(MaxTwin*Fz)
                    gsigs(v, :) = csig(1:round(MaxTwin*Fz))';
                else
                    gsigs(v, 1:length(csig)) = csig(:)';
                end
                ct = (0:length(csig)-1)/Fz;
                plot(ct, csig, 'Color', 0.5*ones(1, 3));
                cEndTid = ctim(2)-ctim(1)+1;
                Moving_end_time_ids(i, k, j, v) = cEndTid;
                plot(ct(cEndTid)*ones(1, 2), [0 60], 'r');
                xlim([0, MaxTwin]);
            end
            Data(i, k, j, :, :) = gsigs;
            assert(sum(~isnan(gsigs(1, :))) > 1);
            cdata = gsigs(:, all(~isnan(gsigs), 1))';
            R2 = corr(cdata+eps*randn(size(cdata))).^2;
            cData(i, k, j) = meanR(R2);
            plot((0:size(gsigs, 2)-1)/Fz, mean(gsigs, 1), 'k');
            if k == 1
                title(sprintf('Speeds:%d(um/s)', IN.Speeds(j)))
                Label_spds{j} = sprintf('%d (um/s)', IN.Speeds(j));
            end
            if j == 1
                if k ==1
                    ylabel(sprintf('Diameters:%d (um)', IN.barWidth(k)));
                else
                    ylabel(sprintf('%d (um)', IN.barWidth(k)));
                end
                Label_dias{k} = sprintf('%d (um)', IN.barWidth(k));
            elseif j == 2
                ylabel('Firing rate (spike/s)')
            end
            ylim([0 180]);
        end
    end
    sgtitle(sprintf('%s Directions:%d', recordingname, IN.direction(i)));
end

%%
figure;
imagesc(log(flipud(squeeze(mean(cData, 1))))); colorbar
yticks(1:length(Label_dias));
yticklabels(fliplr(Label_dias));
xticks(1:length(Label_spds));
xticklabels(Label_spds);
title(sprintf('%s Response quality (log(R2)', recordingname));

%%
dim1_moving_direction = IN.direction;
dim2_bar_width = IN.barWidth;
dim3_speeds = IN.Speeds;
dim4_repeats = 1:IN.nRepeat;
dim5_time = (0:size(Data, 5)-1)/Fz;
saveFileName = sprintf('%s_moving_bar_processed.mat', recordingname);
save(sprintf('./Results/MovingBar/%s', saveFileName), 'dim1_moving_direction', 'dim2_bar_width', 'dim3_speeds',...
    'dim4_repeats', 'dim5_time', 'Data', 'cData', 'Moving_end_time_ids', 'recordingname');