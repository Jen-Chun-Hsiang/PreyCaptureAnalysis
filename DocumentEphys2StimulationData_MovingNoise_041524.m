clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
recording_id = 43;
switch recording_id
    case 1
        recordingname = 'a120723';
        load([ephys_data_folder recordingname '_0010_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
    case 2
        recordingname = 'b120723';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
    case 3
        recordingname = 'c120723';
        load([ephys_data_folder recordingname '_0005_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
    case 4
        recordingname = 'd010224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_2_MovingNoise_1.mat']);
    case 5
        recordingname = 'a030124';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 100;
    case 6
        recordingname = 'c030124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 7
        recordingname = 'd030124';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 8
        recordingname = 'e030124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 9
        recordingname = 'f030124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 150;    
    case 10
        recordingname = 'a030624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 11
        recordingname = 'a030924';
        load([ephys_data_folder recordingname '_0000_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 120;
    case 12
        recordingname = 'b030924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
    case 13
        recordingname = 'c030924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 14
        recordingname = 'd030924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 15
        recordingname = 'a032924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
    case 16
        recordingname = 'b032924';
        load([ephys_data_folder recordingname '_0004_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_004_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 2000;
    case 17
        recordingname = 'c032924';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 18
        recordingname = 'c032924';
        load([ephys_data_folder recordingname '_0003_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
    case 19
        recordingname = 'a040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
    case 25
        recordingname = 'b040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 26
        recordingname = 'c040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 40
        recordingname = 'd040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 41
        recordingname = 'e040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1500;
    case 42
        recordingname = 'f040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1100;
    case 43
        recordingname = 'g040124';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 900;

    case 20
        recordingname = 'a040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 200;
    case 21
        recordingname = 'b040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 22
        recordingname = 'c040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 350;
    case 23
        recordingname = 'd040224';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 24
        recordingname = 'e040224';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 27
        recordingname = 'a040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 150;
    case 28
        recordingname = 'b040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
    case 29
        recordingname = 'c040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 400;
    case 30
        recordingname = 'd040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 31
        recordingname = 'e040524';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 32
        recordingname = 'a040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 600;
    case 33
        recordingname = 'b040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 34
        recordingname = 'c040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 300;
     case 35
        recordingname = 'c040624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 900;
     case 36
        recordingname = 'd040624';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 700;
     case 37
        recordingname = 'e040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 700;
     case 38
        recordingname = 'f040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
     case 39
        recordingname = 'g040624';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 800;
end
%%
OUT = MNS1_OUT;
IN = MNS1_IN;
% OLED.height = 600;
% OLED.width = 800;
% OLED.pixelSize = 2.5; %size of each pixel in micron on retina
clear MNS1_OUT MNS1_IN

assert(abs(length(sectionData)/(OUT.endT-OUT.startT) -10000)<1);
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
%% Use the figure plotted above to determine MinPeakHeight
% perform manual checking and input
%%
sectionData = sectionData-lpass;
%%
Fz = 100;
WinT = [-0.5 0];
% convert trace to spike or firing rate
[pks, locs] = findpeaks(-sectionData, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 40);
figure; hold on
plot(sectionData(1:locs(end)), 'k');
plot(locs(1:end), -pks(1:end), 'rx');
% downsample to Fz
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
%%
UniBlc = unique(OUT.FrmTable(:, 2));
nBlock = length(UniBlc);
bTab = nan(nBlock, 8);
bcount = 1;
for b = 1:length(UniBlc)
    cTab = OUT.FrmTable(OUT.FrmTable(:, 2) == UniBlc(b) , :);
    dT = diff(cTab(:, 1));
    bTab(bcount, :) = [cTab(1, 1), cTab(end, 1), cTab(1, 3),...
        cTab(1, 2), 1/mean(dT),1/min(dT), 1/max(dT),...
        size(cTab, 1)];
    bcount = bcount + 1;
end
bTab(:, 9) = bTab(:, 2)-bTab(:, 1);

%%
NoiBIds = find(bTab(:, 3) == 1);
nRows = ceil(OLED.height/(IN.NoiseGridSize/OLED.pixelSize));          %number of rows of the stimulus array
nCols = ceil(OLED.width/(IN.NoiseGridSize/OLED.pixelSize));
nRows_rs = 600;
nCols_rs = 800;
nFrm = ceil(IN.NoiseFz*IN.NoiseBlockLength);
%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
UniaccuNoiseId = OUT.accuNoiseId;
UniaccuNoiseId(bTab(:, 3) == 11, :) = [];
nProgressSample = 5;
STAmat = zeros(nCols_rs, nRows_rs, ceil(-WinT(1)*Fz));
pSTAmat = zeros(nCols_rs, nRows_rs, ceil(-WinT(1)*Fz));
pstdSTA = zeros(nCols_rs, nRows_rs, nProgressSample);
pBIds = round(linspace(1, length(NoiBIds), nProgressSample+1));
pBIds(1) =  [];
streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
CelData = [];
clear nonlinearindexes
nonlinear_ratio = 0.25;
for i = 1:length(NoiBIds)
    blcFrm = randi(streamU, [0, 1], nCols, nRows, nFrm);
    blcFrm_rs = zeros(800, 600, nFrm);
    for f = 1:nFrm
        cfrm = imresize(blcFrm(:, :, f), [800 600], 'nearest');
        movestep = squeeze(OUT.MovingSteps(NoiBIds(i), f, :));
        xid = movestep(1)+(1:800);
        if movestep(1) >= 0
            xid(xid>800) = [];
            xid_i = 1:length(xid);
        else
            xid_i = find(xid == 1);
            xid_i = xid_i:800;
            xid(xid<1) = [];
        end
        yid = -movestep(2)+(1:600);
        if -movestep(2) >= 0
            yid(yid>600) = [];
            yid_i = 1:length(yid);
        else
            yid_i = find(yid == 1);
            yid_i = yid_i:600;
            yid(yid<1) = [];
        end
        blcFrm_rs(xid, yid, f) = cfrm(xid_i, yid_i);
    end
    lastFrm = squeeze(blcFrm(:, :, end));
    clear blcFrm
    % Verify the noise generator
    verErr = double(UniaccuNoiseId(i, :)>0)-lastFrm(1:size(UniaccuNoiseId, 2));
    if sum(verErr.^2) ~= 0
        error('The white noise reconstruction is wrong!');
    else
        clear lastFrm verErr
        blcFrm_rs = 2*(blcFrm_rs - 0.5);
    end
    
    
    cbTab = bTab(NoiBIds(i), :);
    cData = (find(sig > 0)-1)/Fz;
    %
    fFrm = OUT.FrmTable(:, 1)-OUT.startT;
    %
    frmIds = find(OUT.FrmTable(:, 2) == cbTab(4));
    tFrm = fFrm(frmIds);
    % Spikes (no overlaypping firing rate) and corresponding frame
    nBin = floor((tFrm(end)-tFrm(1))*Fz);
    cnSp = nan(nBin, 4);
    for b = 1:nBin
        cIdx = nBin-(b-1);
        eT = tFrm(end)-(b-1)/Fz;
        sT = eT - 1/Fz;
        cnSp(cIdx, 3) = sT;
        cnSp(cIdx, 4) = eT;
        cnSp(cIdx, 1) = sum(cData>=sT & cData<eT);
        
        dfT = sT - tFrm;
        dfT(dfT<0) = nan;
        [~, cnSp(cIdx, 2)] = nanmin(dfT);
    end
    staCub = zeros(nCols_rs, nRows_rs, ceil(-WinT(1)*Fz));
    talSp = 0;
    bin_time = ceil(-WinT(1)*Fz):nBin;
    nonlinearindexes{i} = randomBinaryVector(bin_time, nonlinear_ratio);
    for b = bin_time(~nonlinearindexes{i})
        if cnSp(b, 1) == 0
            continue
        else
            talSp = talSp + cnSp(b, 1);
            cstaCub = cnSp(b, 1)*blcFrm_rs(:, :, cnSp(ceil(b+WinT(1)*Fz+1):b, 2));
            staCub = staCub + cstaCub;
        end
    end
    staCub = staCub/(talSp+eps);
    STAmat = STAmat + staCub;
    pSTAmat = pSTAmat + staCub;
    if ismember(i, pBIds)
        pBIds_id = find(i == pBIds);
        if pBIds_id == 1
            pstdSTA(:, :, pBIds_id) = std(pSTAmat / i, [], 3);
        else
            pstdSTA(:, :, pBIds_id) = std(pSTAmat / (i-pBIds(pBIds_id-1)), [], 3);
        end
        pSTAmat = zeros(nCols_rs, nRows_rs, ceil(-WinT(1)*Fz));
    end
    CelData.talSp(i) = talSp;
    clear staCub talSp 
end
STAmat = STAmat / length(NoiBIds);
stdSTA = std(STAmat, [], 3);
%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\MEA\Analysis');
close all
%%
smtstdSTA = medfilt2(stdSTA);

[x, y, dist2d] = peak_distance(smtstdSTA);

figure; 
subplot(2, 2, 1);
imagesc(stdSTA');
subplot(2, 2, 2);
imagesc(smtstdSTA');
subplot(2, 2, 3);
imagesc(dist2d');
subplot(2, 2, 4);
scatter(x, y, 5, 'k', 'filled');
%%
figure; 
subplot(1, 2, 1);
a = smtstdSTA';
a = (a-min(a(:)))/range(a(:));
imshow(a);
subplot(1, 2, 2);
imagesc(smtstdSTA'); hold on
plot([400 400], [0 600], '--k');
plot([0 800], [300 300], '--k');
sgtitle(recordingname);
%% get spatial receptive field
minD = 100;
ythr = quantile(y(x>minD), 0.99);
binary_image = smtstdSTA>ythr;
[largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
masked_stdSTA = largest_segment_mask'.*smtstdSTA';
figure; 
subplot(1, 3, 1);
imagesc(binary_image');
subplot(1, 3, 2);
imagesc(largest_segment_mask');
subplot(1, 3, 3);
imagesc(masked_stdSTA);colorbar
%%
smtSTAmat = medfilt3(STAmat);
smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
tRF = reshape(smtSTAmat, [], size(STAmat, 3))'*masked_stdSTA(:);
t = WinT(1):1/Fz:WinT(end);
%%
figure; plot(t(2:end), tRF, 'k');
xlabel('Time (s)');
ylabel('STA contrast');
box off
title(recordingname);

%% Get masked STAmat
csmtSTAmat = permute(smtSTAmat, [2, 1, 3]);
masked_STAmat = nan(size(csmtSTAmat));
for i = 1:size(csmtSTAmat, 3)
    masked_STAmat(:, :, i) = csmtSTAmat(:, :, i).*largest_segment_mask;
end
%%
figure; 
subplot(2, 3, 1);
imagesc(smtstdSTA');
minV = min(pstdSTA(:));
maxV = max(pstdSTA(:));
title('All')
for i = 2:nProgressSample+1
    subplot(2, 3, i)
    imagesc(medfilt2(squeeze(pstdSTA(:, :, i-1)))', [minV, maxV]);
    title(sprintf('Snapshots: %d', i-1));
end
sgtitle(sprintf('%s STA snapshot', recordingname));
%%
% keyboard;
%% Convolute across movies
FRs = [];
PBs = [];
streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
for i = 1:length(NoiBIds)
    blcFrm = randi(streamU, [0, 1], nCols, nRows, nFrm);
    blcFrm_rs = zeros(800, 600, nFrm);
    for f = 1:nFrm
        cfrm = imresize(blcFrm(:, :, f), [800 600], 'nearest');
        movestep = squeeze(OUT.MovingSteps(NoiBIds(i), f, :));
        xid = movestep(1)+(1:800);
        if movestep(1) >= 0
            xid(xid>800) = [];
            xid_i = 1:length(xid);
        else
            xid_i = find(xid == 1);
            xid_i = xid_i:800;
            xid(xid<1) = [];
        end
        yid = -movestep(2)+(1:600);
        if -movestep(2) >= 0
            yid(yid>600) = [];
            yid_i = 1:length(yid);
        else
            yid_i = find(yid == 1);
            yid_i = yid_i:600;
            yid(yid<1) = [];
        end
        blcFrm_rs(xid, yid, f) = cfrm(xid_i, yid_i);
    end
    lastFrm = squeeze(blcFrm(:, :, end));
    clear blcFrm
    % Verify the noise generator
    verErr = double(UniaccuNoiseId(i, :)>0)-lastFrm(1:size(UniaccuNoiseId, 2));
    if sum(verErr.^2) ~= 0
        error('The white noise reconstruction is wrong!');
    else
        clear lastFrm verErr
        blcFrm_rs = 2*(blcFrm_rs - 0.5);
    end
    cbTab = bTab(NoiBIds(i), :);
    cData = (find(sig > 0)-1)/Fz;
    %
    fFrm = OUT.FrmTable(:, 1)-OUT.startT;
    %
    frmIds = find(OUT.FrmTable(:, 2) == cbTab(4));
    tFrm = fFrm(frmIds);
    % Spikes (no overlaypping firing rate) and corresponding frame
    nBin = floor((tFrm(end)-tFrm(1))*Fz);
    cnSp = nan(nBin, 4);
    for b = 1:nBin
        cIdx = nBin-(b-1);
        eT = tFrm(end)-(b-1)/Fz;
        sT = eT - 1/Fz;
        cnSp(cIdx, 3) = sT;
        cnSp(cIdx, 4) = eT;
        cnSp(cIdx, 1) = sum(cData>=sT & cData<eT);
        
        dfT = sT - tFrm;
        dfT(dfT<0) = nan;
        [~, cnSp(cIdx, 2)] = nanmin(dfT);
    end
    bin_time = ceil(-WinT(1)*Fz):nBin;
    bts = bin_time(nonlinearindexes{i}==1);
    FRs = [FRs; cnSp(bts, 1)];
    for b = bts
        cmov = blcFrm_rs(:, :, cnSp(ceil(b+WinT(1)*Fz+1):b, 2));
        PBs = [PBs; mean(cmov.*masked_STAmat, 'all')];
    end
end

figure; scatter(PBs, FRs, 5, 'k', 'filled');
%%
num_data_per_bin = 200; 
[~, sids] = sort(PBs);
num_bin = ceil(length(PBs)/num_data_per_bin);
bin_PBs = nan(num_bin, 1);
bin_FRs = nan(num_bin, 1);
for i = 1:num_bin
    if i < num_bin
        cids = ((i-1)*num_data_per_bin+1):i*num_data_per_bin;
    else
        cids = ((i-1)*num_data_per_bin+1):length(PBs);
    end
    bin_PBs(i) = mean(PBs(sids(cids)));
    bin_FRs(i) = mean(FRs(sids(cids)));
end

figure; scatter(bin_PBs, bin_FRs, 15, 'k', 'filled');
%%
% x = smoothdata(PBs, 'gaussian', 20);
% y = smoothdata(FRs, 'gaussian', 20)*Fz;

x = bin_PBs;
y = bin_FRs*Fz;

x = x./std(x);
% Define the fitting function handle
% 
fit_func = @(p, x) cdf_norm_scaled(x, p(3), p(2), p(1), p(4));
% Initial guesses for scaling parameters
p0 = [30; 1; 1; 0.1];

% Perform the fitting
[pfit, ~] = lsqcurvefit(fit_func, p0, x(:)', y(:)');

% Evaluate the fitted function
y_fit_cdf = cdf_norm_scaled(x, pfit(3), pfit(2), pfit(1), pfit(4));

fit_func = @(p, x) scaledSigmoid(x, p(1), p(2), p(3), p(4));
p0 = [30; -1; 1; 0.1];
% Perform the fitting
[pfit, ~] = lsqcurvefit(fit_func, p0, x(:)', y(:)');
y_fit_sigmoid = scaledSigmoid(x, pfit(1), pfit(2), pfit(3), pfit(4));

figure; hold on
scatter(x, y, 5, 'k', 'filled');
scatter(x, y_fit_cdf, 5, 'b', 'filled');
%scatter(x, y_fit_sigmoid, 5, 'r', 'filled');
xlabel('generator signal');
ylabel('Firing rate');
ylim([0 max(y)]);
title(recordingname);
%%

save([save_data_folder recordingname '.mat'], 'PBs', 'FRs', 'STAmat', 'masked_stdSTA', 'stdSTA',...
    'x', 'y_fit_cdf', 'y_fit_sigmoid', 'masked_STAmat', 'tRF');