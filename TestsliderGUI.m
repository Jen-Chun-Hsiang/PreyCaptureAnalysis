% Create a figure and plot some initial data
clear; close all; clc;
[x1, y1, x2, y2] = sliderGUI()

%%
J = checkerboard(1,2,3);
% Update object positions based on slider values
figure;
hold on;
plot(x1, y1, 'rx'); % Red 'x' for Channel 1
plot(x2, y2, 'go'); % Green 'o' for Channel 2
hold off;
%%
clc; close all;
cgloadlib;
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
IN.mean = 0.5;
IN.NoiseUniqueSeed = 119;
IN.NoiseGridSize = 100; % in um

cgopen(2,0,60,0)
[x1, y1, x2, y2] = WNAlign_sliderGUI(IN, OLED);
fprintf(' x1:%d, y1:%d, \n x2:%d, y2:%d \n', x1, y1, x2, y2);