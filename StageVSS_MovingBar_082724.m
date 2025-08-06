function [OUT] = StageVSS_MovingBar_082624(IN, PROJ, DIO)
%%
speed_n = length(IN.speeds);
width_n = length(IN.Width);
direc_n = length(IN.directions);
bgctr_n = length(IN.BGcontrasts);
speed_ids = 1:speed_n;
width_ids = 1:width_n;
direc_ids = 1:direc_n;
bgctr_ids = 1:bgctr_n;

[X1, X2, X3] = ndgrid(width_ids, speed_ids, direc_ids);
Stims = [X1(:) X2(:) X3(:)];
nStim = size(Stims, 1);
%%
CeterX = IN.Center(1)/PROJ.pixelSize;
CeterY = IN.Center(2)/PROJ.pixelSize;
width_adj = (IN.Width)/PROJ.pixelSize;
height_adj = (IN.Height)/PROJ.pixelSize;
speed_adj = (IN.speeds)/PROJ.pixelSize;
strpoi_adj = IN.startposition/PROJ.pixelSize;
%% color correction
calibration_params = load('NF_IntensityMeasurement_032921_NF_2_002.mat', 'parameters');
calibration_params = calibration_params.parameters;
contrast_bar = ProjIntCalib(calibration_params, IN.spotcontrast);
contrast_bg = ProjIntCalib(calibration_params, IN.BGcontrasts);
%%
time_thread = tic;
if ~PROJ.IsDebug
    pin16high
end
OUT.StartTime = toc(time_thread);
%%
flashlight = stage.builtin.stimuli.Rectangle();
flashlight.size = PROJ.canvasSize;
flashlight.position = PROJ.canvasSize*0.5;
flashlight.color = contrast_bar*ones(1, 3);
presentation_flash = stage.core.Presentation(IN.ISI);
presentation_flash.setBackgroundColor([0, 0, contrast_bg(end)]);
presentation_flash.addStimulus(flashlight);
player_flash = stage.builtin.players.TimeRecordedPrerenderPlayer(presentation_flash, time_thread);
PROJ.client.play(player_flash);
DataTable = nan(IN.repeat*nStim*3, 7); % 6 is the time point
ii = 1;
for b = 1:IN.repeat
    cStims = Stims(randperm(nStim), :);
    cbg = bgctr_ids(randperm(bgctr_n));
    for j = 1:bgctr_n
        for i = 1:nStim
            stim_ids = cStims(i, :);
            fprintf('\t Diameter: %d Speed: %0.2G Direction:%d Contrast:%0.1G \n', IN.Width(stim_ids(1)),...
                IN.speeds(stim_ids(2)), IN.directions(stim_ids(3)), IN.BGcontrasts(cbg(j)));

            % get the start position and moving step
            one_side_distance = strpoi_adj+0.5*width_adj(stim_ids(1));
            [start_x, start_y] = polarToCartesian(one_side_distance,IN.directions(stim_ids(3)));
            start_x = -start_x + PROJ.center_position(1) + CeterX;
            start_y = -start_y + PROJ.center_position(2) + CeterY;
            [mov_x, mov_y] = polarToCartesian(speed_adj(stim_ids(2)), IN.directions(stim_ids(3)));
            stim_duration = 2*one_side_distance/speed_adj(stim_ids(2));

            clc
            fprintf('Moving Spot: %d/%d (repeat) %d/%d \n', b,  IN.repeat, i, nStim);
            % stimulus presenting
            Bar = stage.builtin.stimuli.Rectangle();
            Bar.size  = [width_adj(stim_ids(1)), height_adj];
            Bar.orientation = IN.directions(stim_ids(3));
            Bar.color = [0, 0, contrast_bar];

            % Create a controller to change the aperture's x and y position as a function of time.
            xFunc = @(state)start_x + state.time * mov_x;
            yFunc = @(state)start_y + state.time * mov_y;
            spotPositionController = stage.builtin.controllers.PropertyController(Bar, 'position',...
                @(state)[xFunc(state), yFunc(state)]);

            presentation_movspot = stage.core.Presentation(stim_duration);
            presentation_movspot.setBackgroundColor([0, 0, contrast_bg(cbg(j))]);
            presentation_movspot.addStimulus(Bar);
            presentation_movspot.addController(spotPositionController);
            % Need to add a mask


            % stimulus presenting
            player = stage.builtin.players.TimeRecordedPrerenderPlayer(presentation_movspot, time_thread);
            PROJ.client.play(player);
            play_info = PROJ.client.getPlayInfo();
            DataTable(ii, :) = [b, 1, stim_ids(1), stim_ids(2), stim_ids(3), cbg(j), play_info.flip_times(1)];
            ii = ii +1;
            DataTable(ii, :) = [b, 2, stim_ids(1), stim_ids(2), stim_ids(3), cbg(j), play_info.flip_times(end)];
            ii = ii +1;
            % stimulus offset
            %PROJ.client.play(player_flash);
            %play_info = PROJ.client.getPlayInfo();
            % end of the moving
            epoch_end_time = toc(time_thread);
            while abs(toc(time_thread)-epoch_end_time) < IN.ISI
            end

            DataTable(ii, :) = [b, 0, stim_ids(1), stim_ids(2), stim_ids(3), cbg(j), toc(time_thread)];
            %DataTable(ii, :) = [0, stim_ids(1), stim_ids(2), stim_ids(3), play_info.flip_times(1)];
            ii = ii +1;

        end
    end
end
if ~PROJ.IsDebug
    pin16low
end
OUT.EndTime = toc(time_thread);
OUT.DataTable = DataTable;