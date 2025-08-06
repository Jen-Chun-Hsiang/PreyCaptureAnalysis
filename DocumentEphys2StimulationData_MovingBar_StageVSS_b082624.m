clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar\';

recording_id = 1;
switch recording_id
    case 1
        recordingname = 'a082524';
        load([ephys_data_folder recordingname '_0006_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_006_Retina_1_MovingBar.mat']);
        MinPeakHeight = 500;
    case 2
        recordingname = 'b082524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingBar.mat']);
        MinPeakHeight = 500;
    case 3
        recordingname = 'c082524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingBar.mat']);
        MinPeakHeight = 500;
    case 4
        recordingname = 'd082524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingBar.mat']);
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
num_repeat = IN.repeat;
Speeds = IN.speeds;
barWidth = IN.Width;
directions = IN.directions;
DataTable = OUT.DataTable;
DataTable(:, 7) = round(Fz*(DataTable(:, 6)-OUT.StartTime))+1;

%%
direc_ids = unique(DataTable(:, 4));
speed_ids = unique(DataTable(:, 3));
diame_ids = unique(DataTable(:, 2));
contr_ids = unique(DataTable(:, 5));
Dtab = nan(length(direc_ids)*length(speed_ids)*length(diame_ids)*length(contr_ids), 6);
ii = 1;
for q = 1:length(contr_ids)
    for i = 1:length(direc_ids)
        for j = 1:length(speed_ids)
            for k = 1:length(diame_ids)
                ctime = DataTable(DataTable(:, 4) == direc_ids(i) & DataTable(:, 3) == speed_ids(j) &...
                    DataTable(:, 2) == diame_ids(k) & DataTable(:, 5) == contr_ids(q), 7);
                assert(length(ctime) == 3*num_repeat);
                for v = 1:num_repeat
                    Dtab(ii, :) = [q, i, j, k, ctime(v*3-2), ctime(v*3-1)];
                    ii = ii + 1;
                end
            end
        end
    end
end
Dtab(:, 7) = (Dtab(:, 6)-Dtab(:, 5))/Fz;
%%
MaxTwin = max(Dtab(:, 7))+IN.ISI; % in second (3.2 the max stimulus presentation length)

%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
cData = nan(length(direc_ids), length(contr_ids), length(diame_ids), length(speed_ids));
Data = nan(length(direc_ids), length(contr_ids), length(diame_ids), length(speed_ids), num_repeat, round(MaxTwin*Fz));
Moving_end_time_ids = nan(length(direc_ids), length(contr_ids), length(diame_ids), length(speed_ids),num_repeat);
clear Label_dias Label_spds
close all
for i = 1:length(direc_ids)
    for q = 1:length(contr_ids)
        figure;
        for k = 1:length(diame_ids)
            for j = 1:length(speed_ids)
                subplot(length(diame_ids), length(speed_ids), (k-1)*length(speed_ids)+j); hold on
                cids = find(Dtab(:, 1)==q & Dtab(:, 2)==i & Dtab(:, 4)==k & Dtab(:, 3)==j);
                gsigs = nan(num_repeat, round(MaxTwin*Fz));
                for v = 1:num_repeat
                    ctim = Dtab(cids(v), 5:6);
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
                Data(i, q, k, j, :, :) = gsigs;
                assert(sum(~isnan(gsigs(1, :))) > 1);
                cdata = gsigs(:, all(~isnan(gsigs), 1))';
                R2 = corr(cdata+eps*randn(size(cdata))).^2;
                cData(i, q, k, j) = meanR(R2);
                plot((0:size(gsigs, 2)-1)/Fz, mean(gsigs, 1), 'k');
                if k == 1
                    title(sprintf('Speeds:%d(um/s)', Speeds(j)))
                    Label_spds{j} = sprintf('%d (um/s)', Speeds(j));
                end
                if j == 1
                    if k ==1
                        ylabel(sprintf('Diameters:%d (um)', barWidth(k)));
                    else
                        ylabel(sprintf('%d (um)', barWidth(k)));
                    end
                    Label_dias{k} = sprintf('%d (um)', barWidth(k));
                elseif j == 2
                    ylabel('Firing rate (spike/s)')
                end
                ylim([0 180]);
            end
        end
        sgtitle(sprintf('%s Directions:%d Contrast:%.2G', recordingname, directions(i), IN.BGcontrasts(q)));
    end
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