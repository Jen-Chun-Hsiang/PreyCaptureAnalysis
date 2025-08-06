function OUT = MovingSpot_SciRig_ephys(IN,OLED,DIO)
% time estimation
a = (IN.startposition*2+IN.diameters)./IN.speeds'+IN.ISI;
a = sum(a, 'all')*length(IN.directions)*IN.repeat;

fprintf('estimation time: %0.3G mins \n', a/60);

speed_n = length(IN.speeds);
diameter_n = length(IN.diameters);
direc_n = length(IN.directions);
speed_ids = 1:speed_n;
diameter_ids = 1:diameter_n;
direc_ids = 1:direc_n;

[X, Y, Z] = meshgrid(diameter_ids, speed_ids, direc_ids);
Stims = [X(:) Y(:) Z(:)];
nStim = size(Stims, 1);

rng('shuffle');
%%
CeterX = IN.Center(1)/OLED.pixelSize;
CeterY = IN.Center(2)/OLED.pixelSize;
radius_adj = (IN.diameters*0.5)/OLED.pixelSize;
speed_adj = (IN.speeds)/OLED.pixelSize;
strpoi_adj = IN.startposition/OLED.pixelSize;
%%
calibration_params = load(OLED.ProjIntFile, 'parameters');
calibration_params = calibration_params.parameters;
contrast_spot = ProjIntCalib(calibration_params, IN.spotcontrast);
contrast_bg = ProjIntCalib(calibration_params, IN.bgcontrast);

%%
if ~OLED.IsDebug
    pin16high
end
StartTime = cgflip(contrast_bg,contrast_bg,contrast_bg);
DataTable = nan(IN.repeat*nStim*3, 5); % 3 is the time point
ii = 1;
for b = 1:IN.repeat
    cStims = Stims(randperm(nStim), :);
    for i = 1:nStim
        stim_ids = cStims(i, :);
        % get the start position and moving step
        one_side_distance = strpoi_adj+radius_adj(stim_ids(1));
        [start_x, start_y] = polarToCartesian(one_side_distance,IN.directions(stim_ids(3)));
        start_x = -start_x;
        start_y = -start_y;
        [mov_x, mov_y] = polarToCartesian(speed_adj(stim_ids(3)), IN.directions(stim_ids(3)));
        
        clc
        fprintf('Moving Spot: %d/%d (repeat) %d/%d \n', b,  IN.repeat, i, nStim);
        % stimulus presenting
        Str = 1;
        tic
        while toc*speed_adj(stim_ids(3)) <= one_side_distance*2
            posi_x = start_x+toc*mov_x;
            posi_y = start_y+toc*mov_y;
            cgpencol(contrast_spot,contrast_spot,contrast_spot)
            cgellipse(posi_x,posi_y,radius_adj(stim_ids(1)),radius_adj(stim_ids(1)),'f')
            timestamp = cgflip(contrast_bg,contrast_bg,contrast_bg);
            if Str  == 1
                % start of the moving
                DataTable(ii, :) = [1, stim_ids(1), stim_ids(2), stim_ids(3), timestamp];
                ii = ii +1;
                Str = 0;
                fprintf('\t Diameter: %d Speed: %0.2G Direction:%d \n', IN.diameters(stim_ids(1)),...
                    IN.speeds(stim_ids(2)), IN.directions(stim_ids(3)));
            end
        end
        % end of the moving
        DataTable(ii, :) = [2, stim_ids(1), stim_ids(2), stim_ids(3), timestamp];
        ii = ii +1;
        
        % stimulus offset
        Str = 1;
        tic
        while toc <= IN.ISI
            if Str  == 1
                timestamp = cgflip(contrast_bg,contrast_bg,contrast_bg);
                DataTable(ii, :) = [0, stim_ids(1), stim_ids(2), stim_ids(3), timestamp];
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