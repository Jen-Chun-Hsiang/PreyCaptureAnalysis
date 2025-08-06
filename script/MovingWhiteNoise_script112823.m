%%CLEAR AND INITIALIZE
clear; clc; cgloadlib; 
openparallel; 

%% DEFINE FILENAME UNDER WHICH DATA WILL BE SAVED
Topic = 'Temporal_AlphaRGC';
Day = 'c120723_000_test';
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
OLED.pixelSize = 3.5; %size of each pixel in micron on retina
OLED.A = 0.5857;
OLED.B = 0.4228;
OLED.DRK = 100; %brightness setting in DRK software
OLED.IsDebug = 0;
OLED.ProjIntFile = 'NF_IntensityMeasurement_032921_NF_2_002.mat';

% 1. Moving White Noise
MNS1_IN.uniqueTrigger = 9;
MNS1_IN.mean = 0.5;
MNS1_IN.minT = 10; % 30
MNS1_IN.NoiseUniqueSeed = 119;
MNS1_IN.NoiseRepeatSeed = 78;
MNS1_IN.RepeatProb = 0.125;
MNS1_IN.NoiseGridSize = 100; % in um
MNS1_IN.MovingRange = -50:10:50; % in um
MNS1_IN.NoiseBlockLength = 45; % in second
MNS1_IN.NoiseFz = 10; % 30

% 2. Moving Spot
MS_IN.Center = [0 0]; % (in um)
MS_IN.spotcontrast = 0; % 0: black, 1: white
MS_IN.bgcontrast = 0.8;
MS_IN.speeds = [500 1000 2000 4000 8000]; % (in um)
MS_IN.diameters = [50 100 150 200 300 400 800]; % (in um)
MS_IN.directions = [0 90 180];
MS_IN.startposition = 400; % (in um) how far away from the center
MS_IN.repeat = 2; 
MS_IN.ISI = 0.5;

% 3. Contrast and Stational Spot
SS_IN.Center = [0 0]; % (in um)
SS_IN.diameters = [50 100 150 200 300 400 800]; % (in um)
SS_IN.contrasts = [0:0.1:0.4 0.6:0.1:1]; % (in um)
SS_IN.repeat = 5; 
SS_IN.presentT = 0.5;
SS_IN.ISI = 1;

%% INITIALIZE AND OPEN DISPLAY
if OLED.IsDebug
    cgopen(2,0,60,2)
    OAll.adaptT = 0;
    DIO = nan;
else
    cgopen(2,0,60,2)
end

%% PRESENT MOVING WHITE NOISE
pause(OAll.adaptT)
[MNS1_OUT] = MovingWhiteNoise_SciRig_ephys(MNS1_IN,OLED,DIO);
fileName = sprintf('%s_%s_Retina_%d_MovingNoise_1.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_MovingNoise_1.mat', Topic, Day, RecordId);
end
save([pathName fileName], 'MNS1_IN', 'MNS1_OUT');
clear MNS1_IN MNS1_OUT
fprintf('Moving White Noise 1 ends! %d/%d \n', cSec, OAll.nSec);
cSec = cSec + 1;

%% PRESENT MOVING SPOTS
pause(OAll.adaptT)
[MS_OUT] = MovingSpot_SciRig_ephys(MS_IN,OLED,DIO);
fileName = sprintf('%s_%s_Retina_%d_MovingSpot_1.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_MovingSpot_1.mat', Topic, Day, RecordId);
end
save([pathName fileName], 'MS_IN', 'MS_OUT');
clear MS_IN MS_OUT
fprintf('Moving White Noise 1 ends! %d/%d \n', cSec, OAll.nSec);
cSec = cSec + 1;

%% PRESENT STITIONAL SPOTS
pause(OAll.adaptT)
[SS_OUT] = StationalSpot_SciRig_ephys(SS_IN,OLED,DIO);
fileName = sprintf('%s_%s_Retina_%d_StationalSpot.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_StationalSpot.mat', Topic, Day, RecordId);
end
save([pathName fileName], 'SS_IN', 'SS_OUT');
clear SS_IN SS_OUT
fprintf('Moving White Noise 1 ends! %d/%d \n', cSec, OAll.nSec);
cSec = cSec + 1;
