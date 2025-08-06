clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_name = 'a081024';
load_recording_name = [recording_name '01'];
Fz = 100;
%% get the linear and nonlinear parts
load([load_data_folder load_recording_name '.mat'], 'masked_STAmat', 'PBs', 'FRs');
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_MovingNoise_1.mat'], 'OLED');
% Spatial temporal filter is in "masked_STAmat"

% Nonlinear parts can be derived from "PBs", "FRs"

[nl_fuc, divider]  = getNonlinearFunc(PBs, FRs);
%% Simulation (0) Step full field

x = -1*ones(1, Fz*2);
x(51:150) = 1;
cmov = -1*ones(size(masked_STAmat));
num_time = length(x);
t = (0:num_time-1)/Fz;
resp = nan(1, num_time);
dim1 = size(cmov, 1);
dim2 = size(cmov, 2);
for i = 1:num_time
    resp(i) = mean(cmov.*masked_STAmat, 'all');
    cmov(:, :, 1) = [];
    cmov(:, :, end+1) = x(i)*ones(dim1, dim2);
    clc
    fprintf('Progress ... %d/%d \n', i, num_time);% Da
end
%%
figure; 
subplot(1, 2, 1)
plot(t, resp, 'k');
subplot(1, 2, 2)
plot(t, nl_fuc(resp*0.1/divider), 'k');

%% Simulation (1) varied size spots
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_VariedSizeSpot.mat'], 'SS_IN');
spot_size = SS_IN.diameters;
pix2um = OLED.pixelSize;
[X, Y] = meshgrid(1:size(masked_STAmat, 2), 1:size(masked_STAmat, 1));
masked = sqrt((X-0.5*(size(masked_STAmat, 2)+1)).^2 +  (Y-0.5*(size(masked_STAmat, 1)+1)).^2) < spot_size(6)/pix2um;
%%
figure; 
subplot(1, 2, 1)
imagesc(masked);
subplot(1, 2, 2)
imagesc(masked_STAmat(:, :, 45));

%%
x = zeros(1, Fz*3.5);
x(51:150) = 1;
x(151:300) = -1;
cmov = zeros(size(masked_STAmat));
num_time = length(x);
t = (0:num_time-1)/Fz;
resp = nan(1, num_time);
dim1 = size(cmov, 1);
dim2 = size(cmov, 2);
for i = 1:num_time
    resp(i) = mean(cmov.*masked_STAmat, 'all');
    cmov(:, :, 1) = [];
    cmov(:, :, end+1) = x(i)*masked;
    clc
    fprintf('Progress ... %d/%d \n', i, num_time);
end
%%
figure; 
subplot(1, 2, 1)
plot(t, resp, 'k');
subplot(1, 2, 2)
plot(t, nl_fuc(resp/divider), 'k');
%% Simulation (2) moving bars
BarWidth = 200;
Speed = 2000;
cAng = 0;
a2d = @(x) x*pi/180;
sw = OLED.width;
sh = OLED.height;
diag = sqrt(sw^2+sh^2);
bw = BarWidth/pix2um;
posi = [0 0];
step = Speed.*[cos(a2d(cAng+180)) sin(a2d(cAng+180))];
[X, Y] = meshgrid(1:size(masked_STAmat, 2), 1:size(masked_STAmat, 1));
switch cAng
    case 0
        masked = abs((Y-0.5*(size(masked_STAmat, 1)+1)-posi(1))) < 0.5*bw;
end
%%
figure; 
subplot(1, 2, 1)
imagesc(masked);
subplot(1, 2, 2)
imagesc(masked_STAmat(:, :, 45));

%%
posi = [cos(a2d(cAng))*0.5*diag sin(a2d(cAng))*0.5*diag]+...
            [cos(a2d(cAng))*0.5*bw sin(a2d(cAng))*0.5*bw];
x = ones(1, Fz*3.5);
cmov = -1*ones(size(masked_STAmat));
num_time = length(x);
t = (0:num_time-1)/Fz;
resp = nan(1, num_time);
dim1 = size(cmov, 1);
dim2 = size(cmov, 2);
for i = 1:num_time
    dstep = i*step/Fz;
    move = posi+dstep;
    switch cAng
        case 0
            masked = abs((Y-0.5*(size(masked_STAmat, 1)+1)-move(1))) < 0.5*bw;
    end
    resp(i) = mean(cmov.*masked_STAmat, 'all');
    cmov(:, :, 1) = [];
    cmov(:, :, end+1) = 2*(double(masked)-0.5);
    clc
    fprintf('Progress ... %d/%d \n', i, num_time);
end

%%
figure; 
subplot(1, 2, 1)
plot(t, resp, 'k');
subplot(1, 2, 2)
plot(t, nl_fuc(resp/divider), 'k');
