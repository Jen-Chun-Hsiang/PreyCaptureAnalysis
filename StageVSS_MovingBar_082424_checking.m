clear; clc;
% 
PROJ.IsDebug = 0;

%% DEFINE FILENAME UNDER WHICH DATA WILL BE SAVED
Topic = 'Temporal_AlphaRGC';
Day = 'd082524_001';
RecordId = 1;
cd '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily'
lightcraft_toolbox_dir = 'C:\Users\Kerschensteinerlab\Documents\MATLAB\matlab-lcr';
addpath(genpath(lightcraft_toolbox_dir));
pathName = './Data/Stimulation/';

% OverAll-
OAll.adaptT = 0.1; % 30 in second

% Projector parameters
PROJ.pixelSize = 4.375; %size of each pixel in micron on retina
PROJ.max_led_current = [0 0 50]; % [RGB]
PROJ.target_frame_up_factor = 12;
PROJ.LED_color_str = 'blue'; % 'green'
PROJ.center_position = [415 135]; % in pixel
% get a copy
OAll.pixelSize = PROJ.pixelSize;
OAll.max_led_current = PROJ.max_led_current; % [RGB]
OAll.center_position = PROJ.center_position; % in pixel

 
% 2. Moving Bar (10% contrast)
MB_IN.Center = [0 0]; % (in pixel)
MB_IN.spotcontrast = 1; % 0: black, 1: white
MB_IN.speeds = [500]; % (in um)
MB_IN.Width = [50 100]; % (in um)
MB_IN.BGcontrasts = [0 0.5 0.4 0.3 0.2 0.1]; 
MB_IN.Height = 800; 
MB_IN.directions = [0 90 180];
MB_IN.startposition = 500; % (in um) how far away from the center
MB_IN.repeat = 3; 
MB_IN.ISI = 1;

%%
if ~PROJ.IsDebug
    openparallel;
else
    DIO = 1;
end
%% initiatilize the projector
PROJ.client = stage.core.network.StageClient();
PROJ.client.connect();
PROJ.client.setMonitorGamma(1);

trueCanvasSize = PROJ.client.getCanvasSize();
PROJ.canvasSize = [trueCanvasSize(1) * 2, trueCanvasSize(2)];
PROJ.client.setCanvasProjectionIdentity();
PROJ.client.setCanvasProjectionOrthographic(0, PROJ.canvasSize(1), 0, PROJ.canvasSize(2));

%% Initiate LightCrafter modulation
addpath(genpath('C:\Users\kerschensteinerlab\Documents\MATLAB\matlab-lcr'))
refreshRate = PROJ.client.getMonitorRefreshRate();
PROJ.lightCrafter = LightCrafter4500(refreshRate);
PROJ.lightCrafter.connect();
PROJ.lightCrafter.setMode('pattern');
PROJ.lightCrafter.setLedEnables(true, false, false, false);
 [auto, red, green, blue] = PROJ.lightCrafter.getLedEnables();
 
 %% Set maximum current
[led_currents(1), led_currents(2), led_currents(3)] = PROJ.lightCrafter.getLedCurrents();
% [red_current, green_current, blue_current] = obj.lightCrafter.getLedCurrents();
PROJ.lightCrafter.setLedCurrents(...
    min(led_currents(1), PROJ.max_led_current(1)),...
    min(led_currents(2), PROJ.max_led_current(2)),...
    min(led_currents(3), PROJ.max_led_current(3)));

%% Set up low bit presentation
patternRatesToAttributes = containers.Map('KeyType', 'double', 'ValueType', 'any');

patternRatesToAttributes(1 * refreshRate)  = {8, PROJ.LED_color_str, 1};
patternRatesToAttributes(2 * refreshRate)  = {8, PROJ.LED_color_str, 2};
patternRatesToAttributes(4 * refreshRate)  = {6, PROJ.LED_color_str, 4};
patternRatesToAttributes(6 * refreshRate)  = {4, PROJ.LED_color_str, 6};
patternRatesToAttributes(8 * refreshRate)  = {3, PROJ.LED_color_str, 8};
patternRatesToAttributes(12 * refreshRate) = {2, PROJ.LED_color_str, 12};
patternRatesToAttributes(24 * refreshRate) = {1, PROJ.LED_color_str, 24};

attributes = patternRatesToAttributes(PROJ.target_frame_up_factor*refreshRate);
PROJ.lightCrafter.setPatternAttributes(attributes{:});
            
renderer = stage.builtin.renderers.PatternRenderer(attributes{3}, attributes{1});
PROJ.client.setCanvasRenderer(renderer);

PROJ.current_patternrate = PROJ.lightCrafter.currentPatternRate();


%% PRESENT MOVING SPOTS
pause(OAll.adaptT)
[MB_OUT] = StageVSS_MovingBar_082424(MB_IN,PROJ,DIO);
fileName = sprintf('%s_%s_Retina_%d_MovingBar.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_MovingBar.mat', Topic, Day, RecordId);
end
save([pathName fileName], 'MB_IN', 'MB_OUT', 'OAll');
clear MB_IN MB_OUT
disp('Moving bar ends!');