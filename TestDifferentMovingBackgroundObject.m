seed_id = 100;
exd = 3;
sw = 800;
sh = 400;
SpatialWavelength = 250;
pixelSize = 2.5;
mov_T = 0.3;
stay_T = 0.2;
mov_spd = 2000;
spot_rad = 400; % in um
maxSpd = mov_spd/pixelSize;
maxW = ceil(mov_T*2*maxSpd + sw);
maxH = sh;

stream = RandStream('mrg32k3a','seed',seed_id+1);
sigma = SpatialWavelength/(pixelSize*4);
canvas = imgaussfilt(randn(stream, maxH+round(sigma*exd*2), maxW+round(sigma*exd*2)), sigma);
exdlen = round(exd*sigma);
canvas = canvas((exdlen+1):(size(canvas, 1)-exdlen), (exdlen+1):(size(canvas, 2)-exdlen));
canvas = canvas-min(canvas(:));
canvas = canvas/max(canvas(:));
% canvasW = size(canvas, 1);
% canvasH = size(canvas, 2);
%%
canvas = gaussianblubs(maxW, maxH, SpatialWavelength, pixelSize, seed_id);
figure; imshow(canvas);
hold on
plot([100 200], [100 100], 'y');
plot([100 100], [100 200], 'y');

%%
Stimulus_track_t = repmat([stay_T mov_T], 5, 1)';
Stimulus_speed_bg = [0,  mov_spd;
                     0,  -mov_spd;
                     0, 0;
                     0, -mov_spd;
                     0, mov_spd]';
Stimulus_speed_ob = [0,  0;
                     0,  -mov_spd;
                     0, mov_spd;
                     0, mov_spd;
                     0, -mov_spd]';
%%
Stimulus_track_t = cumsum(Stimulus_track_t(:));
Stimulus_speed_ob = Stimulus_speed_ob(:)/pixelSize;
Stimulus_speed_bg = Stimulus_speed_bg(:)/pixelSize;
spot_rad = spot_rad/pixelSize;
%%
Fz = 100;
Simulation_Time = Stimulus_track_t(end); % in second
simulation_num_frame = round(Simulation_Time*Fz);
extending_length = 100;
gird_axis_w = (-0.5*sw - extending_length):(0.5*sw + extending_length);
gird_axis_h = (-0.5*sh - extending_length):(0.5*sh + extending_length);
[X, Y] = meshgrid(gird_axis_w, gird_axis_h);
cmov = zeros(sh, sw);
spot_center_bg = [0 0];
spot_center_ob = [0 0];
t = 0;
ch = size(canvas, 1);
cw = size(canvas, 2);
close all
figure; 
for i = 1:simulation_num_frame
    % determine the region
    t = t + 1/Fz;
    pid = sum(Stimulus_track_t<t)+1;
    spot_center_bg = [0, spot_center_bg(2) + Stimulus_speed_bg(pid)/Fz];
    spot_center_ob = [0, spot_center_ob(2) + Stimulus_speed_ob(pid)/Fz];
    pids = sqrt((X-spot_center_ob(2)).^2 + (Y-spot_center_ob(1)).^2) < spot_rad;
    x_range = round(ch*0.5-0.5*sh-spot_center_bg(1))+1;
    x_range = x_range:(x_range+sh-1);
    y_range = round(cw*0.5-0.5*sw-spot_center_bg(2))+1;
    y_range = y_range:(y_range+sw-1);
    cmov = canvas(x_range, y_range);
    cmov(pids(extending_length:extending_length+sh-1, extending_length:extending_length+sw-1)) = 0;
    clc
    fprintf('movie %d/%d \n', i, simulation_num_frame);
    spot_center_ob
    spot_center_bg
    imshow(cmov);
    pause(0.2);
end

%%
keyboard;
%% save file
timestamp = datestr(now, 'yyyymmddHHMM');
trajectory_name = sprintf('objectmotion_%s.mat', timestamp);
moviefolder = './Simulation/Trajectories';
save(fullfile(moviefolder, trajectory_name), 'num_dots', 'square_length', 'positions', 'dot_pos');