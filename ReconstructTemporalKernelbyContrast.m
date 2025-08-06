%-----------------------------------------------------------------------
% SCRIPT: computeContrastTemporalKernels.m
%
% This script computes a separate temporal receptive field (tRF) for each
% contrast condition. It does so by:
%  1. Looping over all blocks, grouping them by contrast.
%  2. For each spike in a given block, extracting the preceding
%     num_time_point frames (500 ms = 50 bins at sampling rate Fz),
%     projecting each frame onto a spatial mask (masked_stdSTA),
%     and accumulating those values lagged by time.
%  3. After iterating through all spikes for all blocks of a given
%     contrast, normalizing by the total spike count in that contrast to
%     yield the average temporal kernel (tRF) for that contrast.
%
% REQUIREMENTS / ASSUMPTIONS:
%  - “NoiBIds” is a vector of block indices used above to compute STAmat.
%  - “bTab” is the block table; its 10th column lists the contrast value
%     for each block.
%  - “OUT.contrastPerBlock” gives the contrast for each block index.
%  - “WinT” and “Fz” are as defined in the original STA script:
%       WinT(1) = –0.5 (so that num_time_point = 50 for 500 ms).
%  - A binary (or weighted) spatial mask “masked_stdSTA” of size
%     [nCols_rs × nRows_rs] is already in the workspace.  This mask
%     isolates the RF center; all peripheral pixels should be zero.
%     If you need to threshold “stdSTA” to build this mask, do that
%     beforehand, e.g.:
%         mask = stdSTA;
%         mask(mask < some_threshold) = 0;
%         masked_stdSTA = mask;
%  - The code that reconstructs “blcFrm_rs” per block (as in the STA loop)
%     is essentially identical here, so we repeat it inside each block
%     iteration to rebuild the downsampled stimulus arrays.
%
% OUTPUT:
%  - “contrastVals”: vector of unique contrast values actually used.
%  - “tRFs”: a [num_time_point × numContrasts] matrix; each column is
%     the temporal kernel for that contrast.
%  - “totalSpikesPerContrast”: a [1 × numContrasts] vector, giving the
%     total spike counts that contributed to each tRF.
%-----------------------------------------------------------------------
fit_method = 'least-squares'; % 'least-squares' 'ridge-penalty'
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
% 1) SETUP PARAMETERS
% -------------------
% num_time_point = number of lags = 500 ms before each spike (in bins)
num_time_point = ceil(-WinT(1) * Fz);  
contrastVals = unique( OUT.contrastPerBlock(NoiBIds) );
numContrasts  = numel(contrastVals);

% Preallocate output: each column is the tRF for one contrast value
tRFs = zeros(num_time_point, numContrasts);
% Also store the total spike count per contrast (for normalization / QC)
totalSpikesPerContrast = zeros(1, numContrasts);

% Ensure that the spatial mask is in the workspace:
% masked_stdSTA should be [nCols_rs × nRows_rs].
% If you need to threshold stdSTA first, do that before running this script.
A_cell = cell(1, numContrasts);
B_cell = cell(1, numContrasts);

for c = 1:numContrasts
    A_cell{c} = zeros(num_time_point, num_time_point);
    B_cell{c} = zeros(num_time_point, 1);
end
mask = masked_stdSTA;  
% If masked_stdSTA is already binary or weighted, it can be used directly
% (all peripheral pixels should be zero).

% 2) MAIN LOOP OVER BLOCKS
% ------------------------
for iBlockIdx = 1 : length(NoiBIds)

    clc
    fprintf('TF reconstruct progress... %d/%d \n', iBlockIdx, length(NoiBIds));
    % iBlockIdx is the index in the vector NoiBIds
    blockID = NoiBIds(iBlockIdx);
    thisContrast = OUT.contrastPerBlock(blockID);
    cidx = find( contrastVals == thisContrast );
    
    %% RECONSTRUCT THE STIMULUS FOR THIS BLOCK
    switch lower(IN.intensity_type)
        case 'binary'
            blcFrm  = randi(streamU, [0, 1], nCols, nRows, nFrm) * thisContrast;
        case 'gradient'
            blcFrm  = (rand(streamU, nCols, nRows, nFrm) - 0.5) * thisContrast + 0.5;
        otherwise
            error('Unknown intensity_type in CONTRAST computation.');
    end
    
    % Preallocate resized stimulus arrays
    blcFrm_rs    = zeros(nCols_rs, nRows_rs, nFrm);
    
    for f = 1:nFrm
        cfrm = imresize(blcFrm(:, :, f), [800, 600], 'nearest');
        movestep = squeeze( OUT.MovingSteps(blockID, f, :) );
        xid = movestep(1) + (1:800);
        if movestep(1) >= 0
            xid(xid > 800) = [];
            xid_i = 1:length(xid);
        else
            xid_i = find(xid == 1);
            xid_i = xid_i:800;
            xid(xid < 1) = [];
        end
        
        yid = -movestep(2) + (1:600);
        if -movestep(2) >= 0
            yid(yid > 600) = [];
            yid_i = 1:length(yid);
        else
            yid_i = find(yid == 1);
            yid_i = yid_i:600;
            yid(yid < 1) = [];
        end
        
        blcFrm_rs(xid, yid, f) = cfrm(xid_i, yid_i);
    end
    
    blcFrm_rs    = 2 * (blcFrm_rs    - 0.5);
    
    clear blcFrm cfrm xid yid xid_i yid_i movestep
    
    %% BUILD THE SPIKE-TIME BINS (cnSp) FOR THIS BLOCK
   
    cbTab = bTab(blockID, :);
    cData = (find(sig > 0) - 1) / Fz;   % same as before
    
    fFrm = OUT.FrmTable(:,1) - OUT.startT;
    frmIds = find( OUT.FrmTable(:,2) == cbTab(4) );
    tFrm   = fFrm(frmIds);
    
    nBin = floor(( tFrm(end) - tFrm(1) ) * Fz);
    cnSp = nan(nBin, 4);
    
    for b = 1:nBin
        cIdx = nBin - (b - 1);
        eT = tFrm(end) - (b - 1)/Fz;
        sT = eT - 1/Fz;
        
        cnSp(cIdx, 3) = sT;
        cnSp(cIdx, 4) = eT;
        cnSp(cIdx, 1) = sum( cData >= sT & cData < eT );
        
        dfT = sT - tFrm; 
        dfT(dfT < 0) = nan;
        [~, cnSp(cIdx, 2)] = min(dfT, [], 'omitnan');
    end
    %% LOOP OVER BINS, ACCUMULATE TEMPORAL WEIGHTS
    % -------------------------------------------------------
    % Only bins that start after WinT(1)*Fz (so we have 50 frames before)
    bin_time = (num_time_point : nBin);
    for b = bin_time
        spikeCount = cnSp(b,1);
        if spikeCount == 0
            continue;
        end
    
        %— VECTORIZE the “lag” loop:
        frameIdxs = cnSp(b,2) - (0:(num_time_point-1));  
        valid     = (frameIdxs >= 1);                   
        x_b       = zeros(num_time_point, 1);
    
        if any(valid)
            validFrameIdxs = frameIdxs(valid);                 
            subFrames      = blcFrm_rs(:, :, validFrameIdxs);  % [nRows_rs×nCols_rs×L]
            weighted       = subFrames .* (mask');              % broadcast mask' over 3rd dim
            sums           = squeeze( sum(sum(weighted,1),2) ); % [L×1]
            x_b(valid)     = sums;                             
        end
    
        A_cell{cidx} = A_cell{cidx} + (x_b * x_b');                  % 50×50
        B_cell{cidx} = B_cell{cidx} + (x_b * double(spikeCount));    % 50×1
    
        totalSpikesPerContrast(cidx) = totalSpikesPerContrast(cidx) + spikeCount;
    end

    
    clear blcFrm_rs blcFrm_rs_ds cnSp tFrm fFrm cData
end
save_file_name = [save_recording_name '_contrast_precursorTFs.mat'];
save(fullfile(save_data_folder, save_file_name), 'totalSpikesPerContrast', 'A_cell', 'B_cell');

%% 2) After all blocks: solve for each contrast’s tRF by least‐squares (or ridge)
tRFs = zeros(num_time_point, numContrasts);
lambda = 0.1;   % ridge penalty (set to 0 for pure LS)
for c = 1:numContrasts
    Areg = A_cell{c} + (lambda * eye(num_time_point));
    tRFs(:,c) = Areg \ B_cell{c};
end

%%
save_file_name = [save_recording_name '_contrast_TFs.mat'];
save(fullfile(save_data_folder, save_file_name), 'tRFs', 'lambda');

%% 4) PLOT / SAVE RESULTS
% -----------------------
% Example: Plot all tRFs on one figure
figure; 
hold on;
colors = lines(numContrasts);
timeAxis = (-(num_time_point-1):-1) / Fz * 1000;  % in ms, negative lags
for c = 1:numContrasts
    plot( timeAxis, tRFs(:, c), 'Color', colors(c,:), 'LineWidth', 1.5 );
end
xlabel('Time before spike (ms)');
ylabel('Temporal kernel amplitude');
legend( arrayfun(@num2str, contrastVals, 'UniformOutput', false), ...
        'Location', 'northeast' );
title('Contrast-specific Temporal Receptive Fields (tRF)');
grid on;

% You can also save tRFs to a .mat file if desired:
% save('contrast_tRFs.mat', 'contrastVals', 'tRFs', 'totalSpikesPerContrast');

%-----------------------------------------------------------------------
% END OF SCRIPT
%-----------------------------------------------------------------------
