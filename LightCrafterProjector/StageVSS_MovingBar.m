function [OUT] = StageVSS_MovingBar(IN, PROJ, DIO)
%%
speed_n = length(IN.speeds);
diameter_n = length(IN.diameters);
direc_n = length(IN.directions);
speed_ids = 1:speed_n;
diameter_ids = 1:diameter_n;
direc_ids = 1:direc_n;


[X, Y, Z] = meshgrid(diameter_ids, speed_ids, direc_ids);
Stims = [X(:) Y(:) Z(:)];
nStim = size(Stims, 1);
%%
CeterX = IN.Center(1)/OLED.pixelSize;
CeterY = IN.Center(2)/OLED.pixelSize;
radius_adj = (IN.diameters*0.5)/PROJ.pixelSize;
speed_adj = (IN.speeds)/PROJ.pixelSize;
strpoi_adj = IN.startposition/PROJ.pixelSize;
%%
time_thread = tic;
if ~OLED.IsDebug
    pin16high
end
OUT.StartTime = toc(time_thread);
%%
flashlight = stage.builtin.stimuli.Rectangle();
flashlight.size = PROJ.canvasSize;
flashlight.position = PROJ.canvasSize*0.5;
flashlight.color = 0*ones(1, 3);
presentation_flash = stage.core.Presentation(IN.ISI);
presentation_flash.setBackgroundColor(0);
presentation_flash.addStimulus(flashlight);
player_flash = stage.builtin.players.TimeRecordedPrerenderPlayer(presentation_flash, time_thread);

PROJ.client.play(player_flash);
DataTable = nan(IN.repeat*nStim*3, 5); % 3 is the time point
ii = 1;
for b = 1:IN.repeat
    cStims = Stims(randperm(nStim), :);
    for i = 1:nStim
        stim_ids = cStims(i, :);
        fprintf('\t Diameter: %d Speed: %0.2G Direction:%d \n', IN.diameters(stim_ids(1)),...
            IN.speeds(stim_ids(2)), IN.directions(stim_ids(3)));
       
        % get the start position and moving step
        one_side_distance = strpoi_adj+radius_adj(stim_ids(1));
        [start_x, start_y] = polarToCartesian(one_side_distance,IN.directions(stim_ids(3)));
        start_x = -start_x + PROJ.center_position(1) + CeterX;
        start_y = -start_y + PROJ.center_position(2) + CeterY;
        [mov_x, mov_y] = polarToCartesian(speed_adj(stim_ids(3)), IN.directions(stim_ids(3)));
        stim_duration = 2*one_side_distance/speed_adj(stim_ids(3));
        
        clc
        fprintf('Moving Spot: %d/%d (repeat) %d/%d \n', b,  IN.repeat, i, nStim);
        % stimulus presenting
        spot = stage.builtin.stimuli.Ellipse();
        spot.position = radius_adj(stim_ids(1));
        
        % Create a controller to change the aperture's x and y position as a function of time.
        xFunc = @(state)start_x + state.time * mov_x;
        yFunc = @(state)start_y + state.time * mov_y;
        spotPositionController = stage.builtin.controllers.PropertyController(spot, 'position',...
            @(state)[xFunc(state), yFunc(state)]);
        
        presentation_movspot = stage.core.Presentation(stim_duration);
        presentation_movspot.setBackgroundColor(1);
        presentation_movspot.addStimulus(spot);
        presentation_movspot.addController(spotPositionController);
        % Need to add a mask
        
        
        % stimulus presenting
        player = stage.builtin.players.TimeRecordedPrerenderPlayer(presentation_movspot, time_thread);
        PROJ.client.play(player);
        play_info = PROJ.client.getPlayInfo();
        DataTable(ii, :) = [1, stim_ids(1), stim_ids(2), stim_ids(3), play_info.flip_times(1)];
        DataTable(ii, :) = [2, stim_ids(1), stim_ids(2), stim_ids(3), play_info.flip_times(end)];
        
         % stimulus offset
        PROJ.client.play(player_flash);
        play_info = PROJ.client.getPlayInfo();
        % end of the moving
        DataTable(ii, :) = [0, stim_ids(1), stim_ids(2), stim_ids(3), play_info.flip_times(1)];
        ii = ii +1;
        
    end
end
if ~OLED.IsDebug
    pin16low
end
OUT.EndTime = toc(time_thread);
OUT.DataTable = DataTable;