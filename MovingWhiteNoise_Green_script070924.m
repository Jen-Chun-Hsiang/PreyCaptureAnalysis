 %%CLEAR AND INITIALIZE
clear; clc; cgloadlib;
OLED.IsDebug = 0;
if ~OLED.IsDebug
    openparallel;
end
%% DEFINE FILENAME UNDER WHICH DATA WILL BE SAVED
Topic = 'Temporal_AlphaRGC';
Day = 'temp';
RecordId = 1;
cd '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily'
pathName = './Data/Stimulation/';

%% USER-DEFINED PARAMETERS
% OverAll
OAll.adaptT = 5; % 30 in second
OAll.nSec = 1;
OAll.mean = 0.5;
cSec = 1;

%OLED parameters
OLED.refreshRate = 60;  
OLED.height = 600;
OLED.width = 800;
OLED.pixelSize = 2.5; %size of each pixel in micron on retina
OLED.A = 0.5857;
OLED.B = 0.4228;
OLED.DRK = 100; %brightness setting in DRK software
OLED.ProjIntFile = 'NF_IntensityMeasurement_032921_NF_2_002.mat';

% 1. Moving White Noise
MNS1_IN.uniqueTrigger = 9;
MNS1_IN.mean = 0.5;
MNS1_IN.minT = 12; % 30
MNS1_IN.NoiseUniqueSeed = 119;
MNS1_IN.NoiseRepeatSeed = 78;
MNS1_IN.RepeatProb = 0.125;
MNS1_IN.NoiseGridSize = 100; % in um
MNS1_IN.MovingRange = -50:10:50; % in um
MNS1_IN.NoiseBlockLength = 10; % in second
MNS1_IN.NoiseFz = 10; % 30
MNS1_IN.Color = 'Green'; % Green, UV

% Object Motion 
% ~ 3.x mins
% SNL_IN.uniqueSeed = 101;
% SNL_IN.nRepeat = 2;
% SNL_IN.nTexure = 3;
% SNL_IN.TexureSpatialL = [25 50 100 200 400];
% SNL_IN.nDirection = 1;
% SNL_IN.CycleTime = [
%              0.3 0.7; % stationary
%              0.3 0.7; % Glozbal
%              0.3 0.7; % surround
%              0.3 0.7]; % local
% SNL_IN.Speeds = 1000; % um/s
% SNL_IN.SpatialWavelength = 50; % um
% SNL_IN.SquareLenght = 200; % um
% SNL_IN.RFCenterSize = 300; % um
% SNL_IN.MaskSize = 1000; % um

% 2. Moving Spot
% MS_IN.Center = [0 10]; % (in um)
% MS_IN.spotcontrast = 0; % 0: black, 1: white
% MS_IN.bgcontrast = 1;
% MS_IN.speeds = [500 1000 2000 4000 8000]; % (in um)
% MS_IN.diameters = [50 100 150 200 300 400 800]; % (in um)
% MS_IN.directions = [0 90 180];
% MS_IN.startposition = 400; % (in um) how far away from the center
% MS_IN.repeat = 6; 
% MS_IN.ISI = 0.5;

% 2. Moving Bar
% Current time ~  140
% MB_IN.Center = [0 10]; % (in um)
% MB_IN.stimInt=[1 0 0.6509]; %[backInt stimInt] stimInt= 0 dark bar, 1 light bar
% MB_IN.barWidth = [50 100 150 200 300 400 800]; %width of bars in micron % 50, 100, 300, 500, 800
% MB_IN.speed = [500 1000 2000 4000 8000]; %speed in micron / s (literatures published by Zhou group)
% MB_IN.direction = [0 90 180];
% MB_IN.MaskSize = 2400;
% MB_IN.nRepeats = 6;
% MB_IN.silence = 3; %duration (in (s)) of pause at beginning and end of recording
% MB_IN.ISI = 0.5;

% 3. Varied-size Spot
% SS_IN.Center = [0 10]; % (in um)
% SS_IN.diameters = [50 100 150 200 300 400 800]; % (in um)
% SS_IN.contrasts = 1; % (in um)
% SS_IN.localrepeat = 4; 
% SS_IN.globalrepeat = 2; 
% SS_IN.presentT = [1.2 1.2];
% SS_IN.MaskSize = 1200;
% SS_IN.ISI = 1;

% 3. Varied-contrast Spot
% CS_IN.Center = [0 10]; % (in um)
% CS_IN.diameters = 300; % (in um) % please adjust based on the location Temporal: 300, Nasal: 200/150
% CS_IN.contrasts = 0.2:0.2:1; % (in um)
% CS_IN.localrepeat = 4; 
% CS_IN.globalrepeat = 2; 
% CS_IN.presentT = [1.2 1.2];
% CS_IN.MaskSize = 1200;
% CS_IN.ISI = 1;

%% INITIALIZE AND OPEN DISPLAY
if OLED.IsDebug
    cgopen(2,0,60,3)
    OAll.adaptT = 0;
    DIO = nan;
else
    cgopen(2,0,60,3)
end

%% PRESENT MOVING WHITE NOISE
pause(OAll.adaptT)
[MNS1_OUT] = MovingWhiteNoise_SciRig_ephys_Green(MNS1_IN,OLED,DIO);
fileName = sprintf('%s_%s_Retina_%d_MovingNoise_1.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_MovingNoise_1.mat', Topic, Day, RecordId);
end
save(fullfile(pathName,  fileName), 'MNS1_IN', 'MNS1_OUT', 'OLED');
clear MNS1_IN MNS1_OUT
fprintf('Moving White Noise 1 ends! %d/%d \n', cSec, OAll.nSec);
cSec = cSec + 1;

%% PRESENT OBJECT MOTION
% pause(OAll.adaptT)
% [SNL_OUT] = SpatialNonlinearMotion053124(SNL_IN,OLED,DIO);
% fileName = sprintf('%s_%s_Retina_%d_SpatialNonlinearMotion.mat', Topic, Day, RecordId);
% while exist([pathName fileName], 'file')
%     RecordId = RecordId +1;
%     fileName = sprintf('%s_%s_Retina_%d_SpatialNonlinearMotion.mat', Topic, Day, RecordId);
% end
% save(fullfile(pathName,  fileName), 'SNL_IN', 'SNL_OUT', 'OLED', 'OAll');
% clear SNL_IN SNL_OUT
% fprintf('Spatial nonlinear motion ends! (%d/%d) \n', cSec, OAll.nSec);
% cSec = cSec + 1;

%% PRESENT MOVING SPOTS
% pause(OAll.adaptT)
% [MS_OUT] = MovingSpot_SciRig_ephys(MS_IN,OLED,DIO);
% fileName = sprintf('%s_%s_Retina_%d_MovingSpot_1.mat', Topic, Day, RecordId);
% while exist([pathName fileName], 'file')
%     RecordId = RecordId +1;
%     fileName = sprintf('%s_%s_Retina_%d_MovingSpot_1.mat', Topic, Day, RecordId);
% end
% save([pathName fileName], 'MS_IN', 'MS_OUT', 'OLED');
% clear MS_IN MS_OUT
% fprintf('Moving White Noise 1 ends! %d/%d \n', cSec, OAll.nSec);
% cSec = cSec + 1;

%% PRESENT MOVING SPOTS
% pause(OAll.adaptT)
% [MB_OUT] = MovingBar_SciRig_ephys(MB_IN,OLED,DIO);
% fileName = sprintf('%s_%s_Retina_%d_MovingBar_1.mat', Topic, Day, RecordId);
% while exist([pathName fileName], 'file')
%     RecordId = RecordId +1;
%     fileName = sprintf('%s_%s_Retina_%d_MovingBar_1.mat', Topic, Day, RecordId);
% end
% save([pathName fileName], 'MB_IN', 'MB_OUT', 'OLED');
% clear MB_IN MB_OUT
% fprintf('Moving Bar ends! %d/%d \n', cSec, OAll.nSec);
% cSec = cSec + 1;

%% PRESENT STITIONAL SPOTS
% pause(OAll.adaptT)
% [SS_OUT] = Spot_SciRig_ephys_051724(SS_IN,OLED,DIO);
% fileName = sprintf('%s_%s_Retina_%d_VariedSizeSpot.mat', Topic, Day, RecordId);
% while exist([pathName fileName], 'file')
%     RecordId = RecordId +1;
%     fileName = sprintf('%s_%s_Retina_%d_VariedSizeSpot.mat', Topic, Day, RecordId);
% end
% save(fullfile(pathName,  fileName), 'SS_IN', 'SS_OUT', 'OLED');
% clear SS_IN SS_OUT
% fprintf('Varied-size Spot ends! %d/%d \n', cSec, OAll.nSec);
% cSec = cSec + 1;

%% PRESENT STITIONAL SPOTS
% pause(OAll.adaptT)
% [CS_OUT] = Spot_SciRig_ephys(CS_IN,OLED,DIO);
% fileName = sprintf('%s_%s_Retina_%d_VariedContrastSpot.mat', Topic, Day, RecordId);
% while exist([pathName fileName], 'file')
%     RecordId = RecordId +1;
%     fileName = sprintf('%s_%s_Retina_%d_VariedContrastSpot.mat', Topic, Day, RecordId);
% end
% save([pathName fileName], 'CS_IN', 'CS_OUT', 'OLED');
% clear CS_IN CS_OUT
% fprintf('Varied-contrast Spot ends! %d/%d \n', cSec, OAll.nSec);
% cSec = cSec + 1;
