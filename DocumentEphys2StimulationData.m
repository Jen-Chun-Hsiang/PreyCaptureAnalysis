matFilePath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\EphysStimMatching.mat';
saveEphysData(matFilePath);

%%
EphysFilePath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\';
StructFilePath = ['\\storage1.ris.wustl.' ...
    'edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\EphysStimMatching.mat'];

%% for spliting multiple session (moving white, stationary spot, driftgrating)
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions')
splitEphysData(StructFilePath, EphysFilePath)

%% for 1 section of moving spot
numSection = 1;
splitEphysData(StructFilePath, EphysFilePath, numSection);

%%
OUT = SS_OUT;
assert(abs(length(sectionData)/(OUT.EndTime-OUT.StartTime) -10000)<1);

%%
% convert trace to spike or firing rate
[pks, locs] = findpeaks(sectionData, 'MinPeakHeight', -10, 'MinPeakDistance', 40);
figure; hold on
plot(sectionData(1:locs(1000)), 'k');
plot(locs(1:1000), pks(1:1000), 'rx');
%%
fs = 10000;
signals = zeros(size(sectionData));
signals(locs) = 1;
t = (0:length(sectionData)-1)/fs;
ssigs = smoothdata(signals, 'gaussian', 0.05, 'SamplePoints', t')*fs;
%%
timewindow = 100;
figure; 
subplot(2, 1, 1);
plot(t(1:locs(timewindow)), signals(1:locs(timewindow)), 'k');
subplot(2, 1, 2);
plot(t(1:locs(timewindow)), ssigs(1:locs(timewindow)), 'b');
%%
IN = SS_IN;
dias = IN.diameters;
ndia = length(dias);
ctrs = IN.contrasts;
nctr = length(ctrs);
timewindow = 1.5; % in second

DataTable = OUT.DataTable;
DataTable(:, 5) = round((DataTable(:, 4)-OUT.StartTime)*fs);
assert(DataTable(end, 5)/fs < (OUT.EndTime - OUT.StartTime));

%%
t = 0:(timewindow*fs-1);
%t = -0.1*fs:((timewindow-0.1)*fs-1);
close all
figure; 
for i = 1:ndia
    for j = 1:nctr
        subplot(ndia, nctr, (i-1)*nctr+j);
        cids = find(DataTable(:, 1) == 1 & DataTable(:, 2) == i & DataTable(:, 3) == j);
        fcids = DataTable(cids, 5)' + t';
%         plot(repmat(t', 1, IN.repeat), ssigs(fcids));
        plot(t/fs, mean(ssigs(fcids)', 1), 'k');
        ylim([0 60]);        
        box off
        if j == 1
            ylabel(sprintf('Diameter %d', dias(i)));
        elseif j == 2
            ylabel('Firing rate');
        end
        if i == 1
            title(sprintf('Contrast %.2G', ctrs(j)));
        elseif i == ndia
            xlabel('Time (s)');
        end
    end
end
sgtitle('c120723');