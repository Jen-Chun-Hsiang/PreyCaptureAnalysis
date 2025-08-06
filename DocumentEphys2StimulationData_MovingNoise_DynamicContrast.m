clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
recording_id = 6;
Is_Calculate_STC = 0;
switch recording_id
     case 1
        recordingname = 'a060225';
        save_recording_name = 'a06022501';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_0001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 220;
     case 2
        recordingname = 'a060225';
        save_recording_name = 'a06022502';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_0002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 220;
     case 3
        recordingname = 'c060225';
        save_recording_name = 'c06022501';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_0001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 4
        recordingname = 'c060225';
        save_recording_name = 'c06022502';
        load([ephys_data_folder recordingname '_0002_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_0002_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 5
        recordingname = 'd060225';
        save_recording_name = 'd060225';
        load([ephys_data_folder recordingname '_0004_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_0004_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 1000;
    case 6
        recordingname = 'f060225';
        save_recording_name = 'f060225';
        load([ephys_data_folder recordingname '_0001_1.mat']);
        load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_0001_Retina_1_MovingNoise_1.mat']);
        MinPeakHeight = 500;
end
%%
OUT = MNS1_OUT;
IN = MNS1_IN;
% OLED.height = 600;
% OLED.width = 800;
% OLED.pixelSize = 2.5; %size of each pixel in micron on retina
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
%% Use the figure plotted above to determine MinPeakHeight
% perform manual checking and input
%%
sectionData = sectionData-lpass;
%%
Fz = 100;
WinT = [-0.5 0];
% convert trace to spike or firing rate
[pks, locs] = findpeaks(-sectionData, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 30);
figure; hold on
plot(sectionData(1:locs(end)), 'k');
plot(locs(1:end), -pks(1:end), 'rx');

keyboard;
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
bTab(:, 10) = OUT.contrastPerBlock;
%% ----------------------
% CALIBRATION FILE LOADING
% ----------------------
defaultFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily\Calibration';
C = load(fullfile(defaultFolder,OLED.ProjIntFile));
fprintf('Loaded calibration "%s".\n', fullfile(defaultFolder,OLED.ProjIntFile));

% Single inversion clamp function for both UV and green
inv_clamp = @(y) max(0, min(1, real(C.beta_fit - (1/C.alpha_fit) .* log(C.A_fit ./ (y - C.C_fit + eps) - 1))));



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
num_time_point = ceil(-WinT(1)*Fz);
STAmat = zeros(nCols_rs, nRows_rs, num_time_point);

pSTAmat = zeros(nCols_rs, nRows_rs, num_time_point);
pstdSTA = zeros(nCols_rs, nRows_rs, nProgressSample);
pBIds = round(linspace(1, length(NoiBIds), nProgressSample+1));
pBIds(1) =  [];
streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
CelData = [];
clear nonlinearindexes
nonlinear_ratio = 0.25;
%%
contrastVals = unique( OUT.contrastPerBlock(NoiBIds) );
numContrasts  = numel(contrastVals);
STAmat_ct = zeros(nCols_rs, nRows_rs, num_time_point, numContrasts);
%%
keyboard;
%%
for i = 1:length(NoiBIds)
    clc
    fprintf('STA progress... %d/%d \n', i, length(NoiBIds));
    thisContrast = OUT.contrastPerBlock(NoiBIds(i));
    cidx = find( contrastVals == thisContrast );
    switch lower(IN.intensity_type)
        case 'binary'
            blcFrm  = randi(streamU, [0, 1], nCols, nRows, nFrm)* thisContrast;
        case 'gradient'
            blcFrm  = (rand(streamU, nCols, nRows, nFrm) - 0.5) * thisContrast + 0.5;
        otherwise
            error('The intersity type does not exist.')
    end
    

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
    if i == 1
        lastFrm = inv_clamp(lastFrm);
        % Verify the noise generator
        verErr = double(UniaccuNoiseId(i, :))-lastFrm(1:size(UniaccuNoiseId, 2));
        if sum(verErr.^2) ~= 0
            error('The white noise reconstruction is wrong!');
        end
    end
    clear lastFrm verErr
    blcFrm_rs = 2*(blcFrm_rs - 0.5);
    cbTab = bTab(NoiBIds(i), :);
    cData = (find(sig > 0)-1)/Fz;
    %
    fFrm = OUT.FrmTable(:, 1)-OUT.startT;
    frmIds = find(OUT.FrmTable(:, 2) == cbTab(4));
    tFrm = fFrm(frmIds);
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
    STAmat_ct(:, :, :, cidx) = STAmat_ct(:, :, :, cidx) + staCub;
    
    CelData.talSp(i) = talSp;
    clear staCub talSp
end
STAmat = STAmat / length(NoiBIds);
% STAmat_ct = STAmat_ct / length(NoiBIds);
stdSTA = std(STAmat, [], 3);
%% Get different contrast kernels
contrast_counts = nan(numContrasts, 1);
for i = 1:numContrasts
    contrast_counts(i) = sum(OUT.contrastPerBlock(NoiBIds) == contrastVals(i));
end
%%
STAmat_ct = STAmat_ct./reshape(contrast_counts, 1, 1, 1, []);
%%
colors = lines(numContrasts);
clear contrast_tag
figure; hold on
for i = 1:numContrasts
    contrast_id = i;
    
    smtstdSTA = medfilt2(stdSTA);
    [x, y, dist2d] = peak_distance(smtstdSTA);
    minD = 100;
    ythr = quantile(y(x>minD), 0.99);
    binary_image = smtstdSTA>ythr;
    [largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
    
    smtSTAmat = medfilt3(squeeze(STAmat_ct(:, :, :, contrast_id)));
    masked_stdSTA = largest_segment_mask'.*medfilt2(std(smtSTAmat, [], 3))';
    smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
    tRF = reshape(smtSTAmat, [], size(STAmat, 3))'*masked_stdSTA(:);
    t = WinT(1):1/Fz:WinT(end);
    tRF = (tRF-tRF(end))/range(tRF);
    %
    plot(t(2:end), tRF, 'Color', colors(i, :));
    xlabel('Time (s)');
    ylabel('STA contrast');
    box off
    title(save_recording_name);
    contrast_tag{i} =  sprintf('%2G', contrastVals(i));
end
legend(contrast_tag);

%% Get STC -> use imresize and set covariance size as (10, 10, 50) in the center of (200, 200) pixels
if Is_Calculate_STC
    UniaccuNoiseId = OUT.accuNoiseId;
    UniaccuNoiseId(bTab(:, 3) == 11, :) = [];
    streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
    CelData = [];
    nonlinear_ratio = 0.25;
    all_spikes = 0;
    
    center_region_width_ids = 15:25;
    center_region_height_ids = 10:20;
    num_width_ids = length(center_region_width_ids);
    num_height_ids = length(center_region_height_ids);
    Deviation = nan(num_width_ids, num_height_ids, num_time_point);
    CovMat = zeros(num_width_ids*num_height_ids*num_time_point, num_width_ids*num_height_ids*num_time_point);
    for i = 1:length(NoiBIds)
        thisContrast = OUT.contrastPerBlock(NoiBIds(i));
        blcFrm = randi(streamU, [0, 1], nCols, nRows, nFrm);
        blcFrm_rs = zeros(800, 600, nFrm);
        blcFrm_rs_ds = zeros(downsample_size(1), downsample_size(2), nFrm);
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
            blcFrm_rs_ds(:, :, f) = resize(squeeze(blcFrm_rs(:, :, f)), downsample_size);
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
            [~, cnSp(cIdx, 2)] = min(dfT, [], 'omitnan');
        end
        talSp = 0;
        bin_time = ceil(-WinT(1)*Fz):nBin;
        for b = bin_time(~nonlinearindexes{i})
            talSp = talSp + cnSp(b, 1);
            if cnSp(b, 1) == 0
                continue
            else
                cStimCub = blcFrm_rs_ds(center_region_width_ids, center_region_height_ids, cnSp(ceil(b+WinT(1)*Fz+1):b, 2)) -...
                    STAmat_ds(center_region_width_ids, center_region_height_ids, :);
                for dd = 1:cnSp(b, 1)
                    Deviation = reshape(cStimCub, [], 1)';
                    CovMat = CovMat + Deviation'*Deviation;
                end
            end
        end
        all_spikes = all_spikes + talSp;
    end
    CovMat = CovMat./(all_spikes-1);
    clc
    disp('Finished Covariance Accumulation. Start eigenvalue calculation');

    [eigenvectors, eigenvalues] = eig(CovMat);

    % Convert eigenvalues from matrix to vector
    eigenvalues = diag(eigenvalues);

    % Sort eigenvalues in descending order
    [sorted_eigenvalues, sortIdx] = sort(eigenvalues, 'descend');
    sorted_eigenvectors = eigenvectors(:, sortIdx);

    % Plot the first few eigenvectors (principal components)
    figure;
    for i = 1:3
        subplot(2, 3, i);
        a = reshape(sorted_eigenvectors(:, i), num_width_ids, num_height_ids, num_time_point);
        stdSTA = std(a, [], 3);
        imagesc(stdSTA);
        
        subplot(2, 3, i+3);
        smtSTAmat = permute(a, [2, 1, 3]);
        tRF = reshape(smtSTAmat, [], num_time_point)'*stdSTA(:);
        t = WinT(1):1/Fz:WinT(end);
        plot(t(2:end), tRF, 'k');
        title(sprintf('Principal Component %d', i));
    end
end
%%
% keyboard;
%%
addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\MEA\Analysis');

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
a = (a-min(a(:)))/(max(a(:)-min(a(:))));
imshow(a);
subplot(1, 2, 2);
imagesc(smtstdSTA'); hold on
plot([400 400], [0 600], '--k');
plot([0 800], [300 300], '--k');
sgtitle(save_recording_name);
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
contrast_id = 3;
smtSTAmat = medfilt3(squeeze(STAmat_ct(:, :, :, contrast_id)));
masked_stdSTA = largest_segment_mask'.*medfilt2(std(smtSTAmat, [], 3))';
smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
tRF = reshape(smtSTAmat, [], size(STAmat, 3))'*masked_stdSTA(:);
t = WinT(1):1/Fz:WinT(end);
%
figure; plot(t(2:end), tRF, 'k');
xlabel('Time (s)');
ylabel('STA contrast');
box off
title(save_recording_name);
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
title(save_recording_name);

%% Get masked STAmat
csmtSTAmat = permute(smtSTAmat, [2, 1, 3]);
masked_STAmat = nan(size(csmtSTAmat));
for i = 1:size(csmtSTAmat, 3)
    masked_STAmat(:, :, i) = csmtSTAmat(:, :, i).*largest_segment_mask;
end
%%
% figure;
% subplot(2, 3, 1);
% imagesc(smtstdSTA');
% minV = min(pstdSTA(:));
% maxV = max(pstdSTA(:));
% title('All')
% for i = 2:nProgressSample+1
%     subplot(2, 3, i)
%     imagesc(medfilt2(squeeze(pstdSTA(:, :, i-1)))', [minV, maxV]);
%     title(sprintf('Snapshots: %d', i-1));
% end
% sgtitle(sprintf('%s STA snapshot', save_recording_name));
%%
save([save_data_folder save_recording_name '_STA.mat'], 'STAmat', 'masked_stdSTA', 'stdSTA',...
    'masked_STAmat', 'tRF', 'pstdSTA', 'nProgressSample', 'STAmat_ct');
%%
% keyboard;
%% Convolute across movies
FRs = [];
PBs = [];
streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
for i = 1:length(NoiBIds)
    clc
    fprintf('Contrast progress... %d/%d \n', i, length(NoiBIds));
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
        [~, cnSp(cIdx, 2)] = min(dfT, [], 'omitnan');
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
title(save_recording_name);
%%

save([save_data_folder save_recording_name '.mat'], 'PBs', 'FRs', 'STAmat', 'masked_stdSTA', 'stdSTA',...
    'x', 'y_fit_cdf', 'y_fit_sigmoid', 'masked_STAmat', 'tRF', 'pstdSTA', 'nProgressSample');