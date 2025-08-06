clear; clc;
% 
PROJ.IsDebug = 1;

%% DEFINE FILENAME UNDER WHICH DATA WILL BE SAVED
Topic = 'Temporal_AlphaRGC';
Day = 'test041524_01';
RecordId = 1;
cd '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily'
lightcraft_toolbox_dir = 'C:\Users\Kerschensteinerlab\Documents\MATLAB\matlab-lcr';
addpath(genpath(lightcraft_toolbox_dir));
pathName = './Data/Stimulation/';

% OverAll
OAll.adaptT = 0.1; % 30 in second

% Projector parameters
PROJ.pixelSize = 2.9; %size of each pixel in micron on retina
PROJ.max_led_current = 50;
PROJ.target_frame_up_factor = 12;
PROJ.LED_color_str = 'green';
PROJ.center_position = [415 135]; % in pixel
 
% 2. Moving Bar
MB_IN.Center = [0 0]; % (in pixel)
MB_IN.spotcontrast = 0; % 0: black, 1: white
MB_IN.bgcontrast = 0.8;
MB_IN.speeds = [500 1000 2000 4000 8000]; % (in um)
MB_IN.Width = [50 100 150 200 300 400 800]; % (in um)
MB_IN.Height = 800;
MB_IN.directions = [0 90 180];
MB_IN.startposition = 400; % (in um) how far away from the center
MB_IN.repeat = 5; 
MB_IN.ISI = 0.5;

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
    min(led_currents(1), PROJ.max_led_current),...
    min(led_currents(2), PROJ.max_led_current),...
    min(led_currents(3), PROJ.max_led_current));

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
[MB_OUT] = StageVSS_MovingBar(MB_IN,PROJ,DIO);
fileName = sprintf('%s_%s_Retina_%d_MovingSpot.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_MovingSpot.mat', Topic, Day, RecordId);
end
save([pathName fileName], 'MB_IN', 'MB_OUT');
clear MB_IN MB_OUT
disp('Moving Bar ends!');