function OUT = StationalSpot_SciRig_ephys(IN,OLED,DIO)
%% Time estimation
a = length(IN.diameters)*length(IN.contrasts)*(IN.presentT+IN.ISI);
a = a*IN.repeat;
fprintf('estimation time: %0.3G mins \n', a/60);

contrast_n = length(IN.contrasts);
diameter_n = length(IN.diameters);
contrast_ids = 1:contrast_n;
diameter_ids = 1:diameter_n;

[X, Y] = meshgrid(diameter_ids, contrast_ids);
Stims = [X(:) Y(:)];
nStim = size(Stims, 1);

rng('shuffle');

%%
CeterX = IN.Center(1);
CeterY = IN.Center(2);
radius_adj = (IN.diameters*0.5)/OLED.pixelSize;
%%
calibration_params = load('NF_IntensityMeasurement_032921_NF_2_002.mat', 'parameters');
calibration_params = calibration_params.parameters;
contrast_adj = ProjIntCalib(calibration_params, IN.contrasts);
contrast_mean = ProjIntCalib(calibration_params, 0.5);
%%
if ~OLED.IsDebug
    pin16high
end
StartTime = cgflip(contrast_mean,contrast_mean,contrast_mean);
% wait
tic
while toc <= IN.ISI
end
%
DataTable = nan(IN.repeat*nStim*2, 4);
ii = 1;
for b = 1:IN.repeat
    cStims = Stims(randperm(nStim), :);
    
    for i = 1:nStim
        stim_ids = cStims(i, :);
        clc
        fprintf('Stational Spot: %d/%d (repeat) %d/%d \n', b,  IN.repeat, i, nStim);
        % stimulus presenting
        Str = 1;
        tic
        while toc <= IN.presentT
            if Str  == 1
                cgpencol(contrast_adj(stim_ids(2)),contrast_adj(stim_ids(2)),contrast_adj(stim_ids(2)))
                cgellipse(CeterX,CeterY,radius_adj(stim_ids(1)),radius_adj(stim_ids(1)),'f')
                timestamp = cgflip(contrast_mean,contrast_mean,contrast_mean);
                DataTable(ii, :) = [1, stim_ids(1), stim_ids(2), timestamp];
                ii = ii +1;
                Str = 0;
                fprintf('\t Diameter: %d Contrast: %0.2G \n', IN.diameters(stim_ids(1)),...
                    IN.contrasts(stim_ids(2)));
            end
        end
        % stimulus offset
        Str = 1;
        tic
        while toc <= IN.ISI
            if Str  == 1
                timestamp = cgflip(contrast_mean,contrast_mean,contrast_mean);
                DataTable(ii, :) = [0, stim_ids(1), stim_ids(2), timestamp];
                ii = ii +1;
                Str = 0;
                disp('off set');
            end
        end
    end
end
if ~OLED.IsDebug
    pin16low
end
EndTime = cgflip(contrast_mean,contrast_mean,contrast_mean);
OUT.DataTable = DataTable;
OUT.StartTime = StartTime;
OUT.EndTime = EndTime;

end
