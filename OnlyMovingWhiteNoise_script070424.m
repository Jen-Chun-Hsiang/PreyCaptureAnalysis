 %%CLEAR AND INITIALIZE
clear; clc; cgloadlib;
OLED.IsDebug = 1;
if ~OLED.IsDebug
    openparallel;
end
%% DEFINE FILENAME UNDER WHICH DATA WILL BE SAVED
Topic = 'Temporal_AlphaRGC';
Day = 'temp';
RecordId = 1;
pathName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily\Data\Stimulation\';

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
OLED.pixelSize = 4.35; %size of each pixel in micron on retina
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
MNS1_IN.ColorPos = [0, 0, 0, 0]; %[x1, y1, x2, y2]
MNS1_IN.GreenAttenuation = 0.5;

%% INITIALIZE AND OPEN DISPLAY
if OLED.IsDebug
    cgopen(2,0,60,0, 'Alpha')
    OAll.adaptT = 0;
    DIO = nan;
else
    cgopen(2,0,60,3)
end

%% PRESENT MOVING WHITE NOISE
pause(OAll.adaptT)
[MNS1_OUT] = MovingWhiteNoise_SciRig_ephys_UVGreen(MNS1_IN,OLED,DIO);
fileName = sprintf('%s_%s_Retina_%d_MovingNoise_1.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_MovingNoise_1.mat', Topic, Day, RecordId);
end
save(fullfile(pathName,  fileName), 'MNS1_IN', 'MNS1_OUT', 'OLED');
clear MNS1_IN MNS1_OUT
fprintf('Moving White Noise 1 ends! %d/%d \n', cSec, OAll.nSec);
cSec = cSec + 1;

