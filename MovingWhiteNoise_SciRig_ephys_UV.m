function OUT = MovingWhiteNoise_SciRig_ephys_UV(IN,OLED,DIO)

%
SaveVeriftyNumLength = 20;
% Determine How Many Blocks
nBlock = ceil(IN.minT*60/IN.NoiseBlockLength);

%
rng('shuffle');
RepeatBlockIds = zeros(1, nBlock);
RepeatBlockIds(datasample(1:nBlock, round(nBlock*IN.RepeatProb),'Replace',false)) = 1;
NumRepeat = sum(RepeatBlockIds);
%computing unique stimulus sequence
streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
%computing repeated stimulus sequence
streamR = RandStream('mrg32k3a','seed',IN.NoiseRepeatSeed);

%
nRows = ceil(OLED.height/(IN.NoiseGridSize/OLED.pixelSize));          %number of rows of the stimulus array
nCols = ceil(OLED.width/(IN.NoiseGridSize/OLED.pixelSize));            %number of columns in stimulus array
nFrame = ceil(IN.NoiseFz*IN.NoiseBlockLength);
MovingRange = round(IN.MovingRange/OLED.pixelSize);
MovingSteps = reshape(randsample(MovingRange, nFrame*nBlock*2, 'true'), nBlock, nFrame, []);
MovingSteps_rp = reshape(randsample(MovingRange, nFrame*2, 'true'), nFrame, []);


FrmTable = nan(nBlock*IN.NoiseBlockLength*IN.NoiseFz, 4);
OUT.startT = cgflip;
if ~OLED.IsDebug
    pin16high
end
accuNoiseId = nan(nBlock, SaveVeriftyNumLength);
cFrmId = 1;
for b = 1:nBlock
    clc
    fprintf('White Noise: %d/%d (blocks) \n', b,  nBlock);
    if RepeatBlockIds(b)
        fprintf('\t repeated block: %d/%d \n', sum(RepeatBlockIds(1:b)),  NumRepeat);
        %computing repeated stimulus sequence
        streamR = RandStream('mrg32k3a','seed',IN.NoiseRepeatSeed);
        for i = 1:nFrame
            tic
            Stim_int = randi(streamR, [0, 1], nRows, nCols);
            zero_stim_int = zeros(size(Stim_int));
            %Stim_int = invquadfunc(Stim_int, OLED.A, OLED.B);
            cgloadarray(1, nCols, nRows,[zero_stim_int(:),Stim_int(:),zero_stim_int(:)], OLED.width, OLED.height);
            cgdrawsprite(1,MovingSteps_rp(i, 1),MovingSteps_rp(i, 2))
            while toc < 1/IN.NoiseFz - 0.5*1/60
            end    
            FlipTime = cgflip;            
            FrmTable(cFrmId, :) = [FlipTime, b, 11, i];
            cFrmId = cFrmId + 1;
        end
        accuNoiseId(b, :) = Stim_int(1:SaveVeriftyNumLength);
    else
        for i = 1:nFrame
            tic
            Stim_int = randi(streamU, [0, 1], nRows, nCols);
            zero_stim_int = zeros(size(Stim_int));
            %Stim_int = invquadfunc(Stim_int, OLED.A, OLED.B);
            cgloadarray(1, nCols, nRows,[zero_stim_int(:),Stim_int(:),zero_stim_int(:)], OLED.width, OLED.height);
            cgdrawsprite(1,MovingSteps(b, i, 1),MovingSteps(b, i, 2))
            while toc < 1/IN.NoiseFz - 0.5*1/60
            end   
            FlipTime = cgflip;            
            FrmTable(cFrmId, :) = [FlipTime, b, 1, i];
            cFrmId = cFrmId + 1;
        end
        accuNoiseId(b, :) = Stim_int(1:SaveVeriftyNumLength);
    end
end
OUT.endT = cgflip;
pin16low
cgflip(0,0,0)
cgflip
OUT.FrmTable = FrmTable;
OUT.accuNoiseId = accuNoiseId;
OUT.nBlock = nBlock;
OUT.RepeatBlockIds = RepeatBlockIds;
OUT.MovingSteps = MovingSteps;
OUT.MovingSteps_rp = MovingSteps_rp;

end