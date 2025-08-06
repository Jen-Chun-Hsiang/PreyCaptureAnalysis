% The goal of this script is to understand how a neuron is trigger by
% stimulus. The neuron composed of a multi-lobed receptive field
clear; close all; clc
%% Parameters
% a040124, b040224
file_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
save_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\SingleGaussianSplit';
file_name = 'a040124';
load(fullfile(file_folder, [file_name, '.mat']));

%% Load STA
InitThresQuant = 0.98;
num_gauss = 6;
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

%%
keyboard;
%%

image = stdSTA';
initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];
initial_params = repmat(initial_params, num_gauss, 1);
initial_params(:, 1:2) = initial_params(:, 1:2) + 10*rand(num_gauss, 2);
initial_params = initial_params';

%objective_function = @(params) gaussian_difference(params, image, num_gauss);
objective_function = @(params) 1-corr(image(:), reshape(gaussian_multi(params, image, num_gauss), [], 1));
% objective_function = @(params) sum((image(:) - reshape(gaussian_multi(params, image, num_gauss), [], 1)).^2);
options.MaxFunEvals = 600*length(initial_params(:));
[optimal_params,fval,exitflag,output]= fminsearch(objective_function, initial_params);

gaussian_model = gaussian_multi(optimal_params, image, num_gauss);
subplot(1, 2, 2);
imagesc(gaussian_model');
hold on
plot(300, 400, 'rx');
% scatter(optimal_params(2, :), optimal_params(1, :), 25, 'r', 'filled');
for i = 1:size(optimal_params, 2)
    text(optimal_params(2, i), optimal_params(1, i), sprintf('%d', i));
end
cost = objective_function(optimal_params);

save_file_name = sprintf('%s_dot%d.mat', file_name, num_gauss);
save(fullfile(save_folder, save_file_name), 'num_gauss', 'image', 'optimal_params', 'cost');
%% Separate into different individual lobed receptive fields
% Removed the lobe that is outside the mask
removeids = zeros(num_gauss, 1);
for i = 1:num_gauss
    x = round(optimal_params(2, i));
    y = round(optimal_params(1, i));
    if x < 0 || x > size(masked_stdSTA, 1)
        removeids(i) = 1;
        continue
    end
    if y < 0 || y > size(masked_stdSTA, 2)
        removeids(i) = 1;
        continue
    end
    if masked_stdSTA(x, y) == 0
        removeids(i) = 1;
    end
end
optimal_params = optimal_params(:, removeids==0);
num_gauss = size(optimal_params, 2);
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
    c = gaussian2d(X, Y, optimal_params(:, i))';
    lobed_STAs(:, :, :, i) = repmat(c, 1, 1, length(tRF), 1).*reshape(tRF, 1, 1, length(tRF));
end
% Get temporal filter
%%
[X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
figure;
for i = 1:num_gauss
    subplot(1, num_gauss, i);
    c = gaussian2d(X, Y, optimal_params(:, i))';
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
MinPeakHeight = 300;

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
Conv_blobSTAs = [];
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
    keyboard;
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
    staCub = zeros(nCols_rs, nRows_rs, ceil(-WinT(1)*Fz));
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
    Conv_blobSTAs = cat(1, Conv_blobSTAs, temp_Conv_blobSTAs);
end
save_file_name = sprintf('%s_dot%d_response_profile.mat', file_name, num_gauss);
save(fullfile(save_folder, save_file_name), 'num_gauss', 'Conv_blobSTAs', 'optimal_params', 'Select_Spikes', 'Fz');
%% Examine the correlation between spikes and individual lobed receptive field
figure; hold on
num_bin = 50;
Colors = parula(num_gauss);
for i = 1:num_gauss
    % subplot(2, 2, i);
    a = Conv_blobSTAs(:, i);
    xbins = linspace(min(a(:)), max(a(:)), num_bin);
    x = nan(num_bin-1, 2);
    for j = 1:num_bin-1
        cids = a>= xbins(j) & a<xbins(j+1);
        x(j, 1) = mean(a(cids));
        x(j, 2) = mean(Select_Spikes(cids));
    end
    plot(x(:, 1), x(:, 2)*Fz, 'Color', Colors(i, :));
    % plot(x(:, 1), x(:, 2)*Fz, 'k'); % 'Color', Colors(i, :));
    xlabel('Blob intensity');
    ylabel('Firing rate (spike/s)');
    % ylim([0 100]);
    % xlim([-2e7 2e7])
    % title(sprintf('Lobe %d', i));
    box off
end

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
            cids = cids & Conv_blobSTAs(:, j) <0;
        else
            cids = cids & Conv_blobSTAs(:, j) >=0;
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
yticks(1:3:16);
yticklabels({'16', '13', '10', '7', '4', '1'});
box off
xlabel('Firing rate (spikes/s)');
subplot(1, 4, 2);
imagesc(D_comb);
colorbar
title('Firing rate map');
subplot(1, 4, 4);
barh(flipud(D_count));
ylim([0.5 length(D_count)+0.5]);
yticks(1:3:16);
yticklabels({'16', '13', '10', '7', '4', '1'});
box off
xlabel('Time bin counts');

