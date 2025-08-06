function OUT = Spot_SciRig_ephys_051724(IN,OLED,DIO)
%% Time estimation
a = length(IN.diameters)*length(IN.contrasts)*sum(IN.presentT)+IN.ISI;
a = a*IN.localrepeat*IN.globalrepeat;
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
radius_adj = IN.diameters/OLED.pixelSize;
maskdia_adj = IN.MaskSize/OLED.pixelSize;
%%
calibration_params = load('NF_IntensityMeasurement_032921_NF_2_002.mat', 'parameters');
calibration_params = calibration_params.parameters;
contrast_adj_ON = ProjIntCalib(calibration_params, 0.5+(IN.contrasts/2));
contrast_adj_OFF = ProjIntCalib(calibration_params, 0.5-(IN.contrasts/2));
contrast_mean = ProjIntCalib(calibration_params, 0.5);
%%
if ~OLED.IsDebug
    pin16high
end
cgpencol(contrast_mean, contrast_mean, contrast_mean)
cgellipse(CeterX,CeterY,maskdia_adj,maskdia_adj,'f')
StartTime = cgflip(0,0,0);
% wait
tic
while toc <= IN.ISI
end
%
DataTable = nan(IN.localrepeat*IN.globalrepeat*nStim*2, 7);
ii = 1;
for b = 1:IN.globalrepeat
    cStims = Stims(randperm(nStim), :);
    for i = 1:nStim
        stim_ids = cStims(i, :);
        for r = 1:IN.localrepeat
            clc
            fprintf('Repeats: %d/%d, %d/%d, Stational Spot: %d/%d \n',...
                b,  IN.globalrepeat, r, IN.localrepeat, i, nStim);
            
            for q = 1:2
                switch mod(b+q, 2)
                    case 0
                        % ON stimulus presenting
                        Str = 1;
                        tic
                        while toc <= IN.presentT(1)
                            if Str  == 1
                                cgpencol(contrast_mean, contrast_mean, contrast_mean)
                                cgellipse(CeterX,CeterY,maskdia_adj,maskdia_adj,'f')
                                cgpencol(contrast_adj_ON(stim_ids(2)),contrast_adj_ON(stim_ids(2)),contrast_adj_ON(stim_ids(2)))
                                cgellipse(CeterX,CeterY,radius_adj(stim_ids(1)),radius_adj(stim_ids(1)),'f')
                                timestamp = cgflip(0,0,0);
                                DataTable(ii, :) = [b, i, r, 1, stim_ids(1), stim_ids(2), timestamp];
                                ii = ii +1;
                                Str = 0;
                                fprintf('\t Diameter: %d Contrast: %0.2G \n', IN.diameters(stim_ids(1)),...
                                    IN.contrasts(stim_ids(2)));
                            end
                        end
                    case 1
                        % OFF stimulus presenting
                        Str = 1;
                        tic
                        while toc <= IN.presentT(2)
                            if Str  == 1
                                cgpencol(contrast_mean, contrast_mean, contrast_mean)
                                cgellipse(CeterX,CeterY,maskdia_adj,maskdia_adj,'f')
                                cgpencol(contrast_adj_OFF(stim_ids(2)),contrast_adj_OFF(stim_ids(2)),contrast_adj_OFF(stim_ids(2)))
                                cgellipse(CeterX,CeterY,radius_adj(stim_ids(1)),radius_adj(stim_ids(1)),'f')
                                timestamp = cgflip(0,0,0);
                                DataTable(ii, :) = [b, i, r, 0, stim_ids(1), stim_ids(2), timestamp];
                                ii = ii +1;
                                Str = 0;
                            end
                        end
                end
            end
        end
        cgpencol(contrast_mean, contrast_mean, contrast_mean)
        cgellipse(CeterX,CeterY,maskdia_adj,maskdia_adj,'f')
        cgflip(0,0,0);
        disp('ISI');
        tic
        while toc <= IN.ISI
        end
    end
end

if ~OLED.IsDebug
    pin16low
end
cgpencol(contrast_mean, contrast_mean, contrast_mean)
cgellipse(CeterX,CeterY,maskdia_adj,maskdia_adj,'f')
EndTime = cgflip(0,0,0);
OUT.DataTable = DataTable;
OUT.StartTime = StartTime;
OUT.EndTime = EndTime;

end
