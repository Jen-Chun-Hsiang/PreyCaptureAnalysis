%%CLEAR AND INITIALIZE
clear; clc; cgloadlib; 

%% DEFINE FILENAME UNDER WHICH DATA WILL BE SAVED
Topic = 'MEA_Dunnart';
Day = '20230216';
RecordId = 1;
% cd 'C:\Users\kerschensteinerlab\Documents\MATLAB\Scripts\VisualStimulation\currentScripts\stimScript_JH'
pathName = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily\Data';

%% USER-DEFINED PARAMETERS
%OLED parameters
OLED.refreshRate = 60;  
OLED.height = 600;
OLED.width = 800;
OLED.pixelSize = 2.5; %size of each pixel in micron on retina
OLED.A = 0.5766;
OLED.B = 0.4297 ;
OLED.DRK = 98; %brightness setting in DRK software
OLED.IsDebug = 1;

% OverAll
OAll.adaptT = 5; % 30 in second
OAll.spontT = 300; % 300 in second
OAll.nSec = 10;
OAll.mean = 0.5;
cSec = 1;

% Object Motion 
% ~ 21 mins
SNL_IN.uniqueTrigger = 10;
SNL_IN.uniqueSeed = 101;
SNL_IN.nRepeat = 2;
SNL_IN.nTexure = 3;
SNL_IN.TexureSpatialL = [25 50 100 200 400];
SNL_IN.nDirection = 1;
SNL_IN.CycleTime = [
             0.3 0.2; % stationary
             0.3 0.7; % Glozbal
             0.3 0.7; % surround
             0.3 0.7]; % local
SNL_IN.Speeds = 500; % um/s
SNL_IN.SpatialWavelength = 50; % um
SNL_IN.SquareLenght = 200; % um
SNL_IN.RFCenterSize = 300; % um
SNL_IN.MaskSize = 1000; % um

%% INITIALIZE AND OPEN DISPLAY
disp(['switch on OLED and set brightness to ' num2str(OLED.DRK(1)) ' in software'])

if OLED.IsDebug
    cgopen(2,0,60,0)
    OAll.adaptT = 0;
    DIO = 0;
else
    openparallel; 
    cgopen(2,0,60,2)
    displayCorrect(OLED, OAll.mean)
%     OAll.spontT = 0;
%     OAll.adaptT = 0;
end

%% PRESENT OBJECT MOTION
pause(OAll.adaptT)
[SNL_OUT] = SpatialNonlinearMotion053124(SNL_IN,OLED,DIO);
fileName = sprintf('%s_%s_Retina_%d_DiffMotionGrid_1.mat', Topic, Day, RecordId);
while exist([pathName fileName], 'file')
    RecordId = RecordId +1;
    fileName = sprintf('%s_%s_Retina_%d_DiffMotionGrid_1.mat', Topic, Day, RecordId);
end
save([pathName fileName], 'SNL_IN', 'SNL_OUT', 'OLED', 'OAll');
clear SNL_IN SNL_OUT
fprintf('Spatial nonlinear ends! (%d/%d) \n', cSec, OAll.nSec);
cSec = cSec + 1;