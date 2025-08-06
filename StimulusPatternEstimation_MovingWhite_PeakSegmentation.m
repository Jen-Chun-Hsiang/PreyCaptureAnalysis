% The goal of this script is to understand how a neuron is trigger by
% stimulus. The neuron composed of a multi-lobed receptive field
clear; close all; clc
%% Parameters
% a040124, b040224
file_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\SingleGaussianSplit';
file_name = 'a040624';
load(fullfile(file_folder, [file_name, '.mat']));

%% Load STA
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
InitThresQuant = 0.98;
num_gauss = 7;
close all; clc
figure;
subplot(1, 3, 1);
imagesc(stdSTA);colorbar
hold on
plot(300, 400, 'rx');

img = stdSTA;
Threshold = quantile(img(:), InitThresQuant);
if Threshold == 0
    Threshold = mean(img(:));
end
% Initialize the potential synapses
Cent = FastPeakFind(img, Threshold);
Cent = reshape(Cent, 2, [])';
index = sub2ind(size(img), Cent(:, 2), Cent(:, 1));
rmids = img(index) < Threshold;
Cent(rmids, :) = [];

hold on; scatter(Cent(:, 1), Cent(:, 2), 'wo');
%%
subplot(1, 3, 2);
imagesc(stdSTA);colorbar
hold on
plot(300, 400, 'rx');
index = sub2ind(size(img), Cent(:, 2), Cent(:, 1));
Cent = removeNeighborPeak(Cent, img(index), num_gauss);

hold on; scatter(Cent(:, 1), Cent(:, 2), 'wo');
%%
MaxY = size(img, 1);
MaxX = size(img, 2);
[X, Y] = meshgrid(1:MaxX, 1:MaxY);
NumCent = size(Cent, 1);
DistCent = nan(numel(img), NumCent);
for i = 1:NumCent
    DistCent(:, i) = sqrt((X(:)-Cent(i, 1)).^2+(Y(:)-Cent(i, 2)).^2);
end
[~, Idx] = min(DistCent, [], 2);
ClusterMap = reshape(Idx, size(img, 1), []);
ClusterMap(img < Threshold*0.65) = 0;
MapMask = zeros(size(ClusterMap));
for i = 1:NumCent
    CC = bwconncomp(ClusterMap == i, 4);
    CentComp = cellfun(@(x) ismember(sub2ind(size(img), Cent(i, 2), Cent(i, 1)), x), CC.PixelIdxList);
    if any(CentComp)
        MapMask(CC.PixelIdxList{CentComp == 1}) = i;
    end
end
clear CC Idx
subplot(1, 3, 3);
imagesc(MapMask);colorbar;


%% Separate into different individual lobed receptive fields
%%
Fz = 100;
WinT = [-0.5 0];
smtSTAmat = medfilt3(STAmat);
smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
tRF = reshape(smtSTAmat, [], size(STAmat, 3))'*masked_stdSTA(:);
t = WinT(1):1/Fz:WinT(end);
%
figure; plot(t(2:end), tRF, 'k');
xlabel('Time (s)');
ylabel('STA contrast');
box off
%% Visualize and reconstruct the STA
[X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
dims = size(STAmat);
lobed_STAs = nan(dims(1), dims(2), dims(3), num_gauss);
for i = 1:num_gauss
    c = double(MapMask == i).*stdSTA;
    lobed_STAs(:, :, :, i) = repmat(c, 1, 1, length(tRF), 1).*reshape(tRF, 1, 1, length(tRF));
end
% Get temporal filter
%%
[X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
figure;
for i = 1:num_gauss
    subplot(1, num_gauss, i);
    c = double(MapMask == i).*stdSTA;
    imagesc(c);
    hold on
    plot(300, 400, 'rx');
    title(sprintf('Lobe %d', i));
    box off
end
%% Convoluted with each frame can calculate that how each lobed receptive field is related to the spikes
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recordingname = file_name;
load([ephys_data_folder recordingname '_0001_1.mat']);
load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingNoise_1.mat']);
MinPeakHeight = 1000;

%%
OUT = MNS1_OUT;
IN = MNS1_IN;
clear MNS1_OUT MNS1_IN

assert(abs(length(sectionData)/(OUT.endT-OUT.startT) -10000)<1);
%%
if ~exist('save_recording_name', 'var')
    save_recording_name = recordingname;
end
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
nFrm = ceil(IN.NoiseFz*IN.NoiseBlockLength);

%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
UniaccuNoiseId = OUT.accuNoiseId;
UniaccuNoiseId(bTab(:, 3) == 11, :) = [];
num_time_point = ceil(-WinT(1)*Fz);
streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
clear nonlinearindexes
nonlinear_ratio = 0.01;

%%
Conv_lobeSTAs = [];
Select_Spikes = [];
for i = 1:length(NoiBIds)
    clc
    fprintf('Conv progress... %d/%d \n', i, length(NoiBIds));
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
        [~, cnSp(cIdx, 2)] = min(dfT, [], 'omitnan'); % nanmin(dfT)
    end
    
    talSp = 0;
    bin_time = ceil(-WinT(1)*Fz):nBin;
    nonlinearindexes{i} = randomBinaryVector(bin_time, nonlinear_ratio);
    bts = bin_time(~nonlinearindexes{i});
    temp_Conv_blobSTAs = nan(length(bts), num_gauss);
    Select_Spikes = [Select_Spikes; cnSp(bts, 1)];
    for b = 1:length(bts)
        cmov = blcFrm_rs(:, :, cnSp(ceil(bts(b)+WinT(1)*Fz+1):bts(b), 2));
        if b == 1 % make sure the size is matched
            dims = size(lobed_STAs);
            assert(sum(size(cmov)-dims(1:3))== 0);
        end

        for k = 1:num_gauss
            temp_Conv_blobSTAs(b, k) = cmov(:)'*reshape(lobed_STAs(:, :, :, k), [], 1);
        end
    end
    Conv_lobeSTAs = cat(1, Conv_lobeSTAs, temp_Conv_blobSTAs);
end
save_file_name = sprintf('%s_dot%d__peak_response_profile.mat', file_name, num_gauss);
save(fullfile(save_folder, save_file_name), 'num_gauss', 'Conv_lobeSTAs', 'MapMask', 'Select_Spikes', 'Fz');
%%
%% Examine the correlation between spikes and individual lobed receptive field

figure; hold on
num_bin = 15;
Colors = parula(num_gauss);
for i = 1:num_gauss
    % subplot(2, 2, i);
    a = Conv_lobeSTAs(:, i);
    % xbins = linspace(min(a(:)), max(a(:)), num_bin);
    xbins = equalspacingbin(a, num_bin);
    x = nan(num_bin-1, 2);
    for j = 1:num_bin
        cids = a>= xbins(j, 1) & a<xbins(j, 2);
        x(j, 1) = mean(a(cids));
        x(j, 2) = mean(Select_Spikes(cids));
    end
    plot(x(:, 1), x(:, 2)*Fz, 'Color', Colors(i, :));
    % plot(x(:, 1), x(:, 2)*Fz, 'k'); % 'Color', Colors(i, :));
    xlabel('Lobe intensity');
    ylabel('Firing rate (spike/s)');
    % ylim([0 100]);
    % xlim([-2e7 2e7])
    % title(sprintf('Lobe %d', i));
    box off
end
%%
xbins = linspace(-1e5, 1e5, 50);
Colors = parula(num_gauss);
figure; hold on
for i = 1:num_gauss
    h1 = histogram(Conv_lobeSTAs(:, i), xbins);
    h1.FaceColor = Colors(i, :);
    h1.Normalization = 'Probability';
    h1.EdgeColor = 'w';
    h1.FaceAlpha = 0.5;
end
xlabel('Lobe intensity');
ylabel('Probability');
%% Examine how the combination of lobed receptive field patterns contributed to the spikes
% (normalize to the trials of each pattern)
combinations = dec2bin(0:2^num_gauss-1) - '0';
num_comb = size(combinations, 1);
D_comb = nan(num_comb, 1);
D_count = nan(num_comb, 1);
for i = 1:num_comb
    cids = ~isnan(Select_Spikes);
    for j = 1:num_gauss
        if combinations(i, j) == 0
            cids = cids & Conv_lobeSTAs(:, j) <0;
        else
            cids = cids & Conv_lobeSTAs(:, j) >=0;
        end
    end
    D_comb(i) = mean(Select_Spikes(cids));
    D_count(i) = sum(cids);
end
%%
figure;
subplot(1, 4, 1);
imagesc(combinations);
box off
title('Patterns');
subplot(1, 4, 3);
barh(flipud(D_comb));
ylim([0.5 length(D_comb)+0.5]);
yticks([]);
% yticklabels({'60', '50', '40', '30', '20', '10'});
box off
xlabel('Firing rate (spikes/s)');
subplot(1, 4, 2);
imagesc(D_comb);
colorbar
title('Firing rate map');
subplot(1, 4, 4);
barh(flipud(D_count));
ylim([0.5 length(D_count)+0.5]);
yticks([]);
% yticklabels({'16', '13', '10', '7', '4', '1'});
box off
xlabel('Time bin counts');

%%
X = combinations;
y = D_comb;
s = sum(combinations, 2);
unique_s = unique(s);
num_s = length(unique_s);
[weights, y_pred] = calculateWeights(X, y);
Colors = parula(num_s);
figure; 
subplot(1, 2, 1); hold on
for i = 1:num_s
    cids = s == unique_s(i);
    scatter(y_pred(cids)*Fz, y(cids)*Fz, 15, Colors(i, :), 'filled');
end
plot([0.1 0.7]*Fz, [0.1 0.7]*Fz, '--k');
ylabel('Observed firing rate (spike/s)');
xlabel('Predicted firing rate (spike/s)');
xlim([10 70]);
ylim([10 70]);
subplot(1, 2, 2); hold on
imagesc(s); colorbar
%%
Colors = parula(num_gauss);
figure; hold on
for i = 1:num_gauss
    bh = bar(i, weights(1+i));
    bh.FaceColor = Colors(i, :);
    bh.EdgeColor = 'w';
end
ylabel('Weights');
box off