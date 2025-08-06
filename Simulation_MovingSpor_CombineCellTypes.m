clear; close all; clc
%%
OFFCellFileName = 'OFFRGCs_Connect1_042824_trajectory05.mat';
Compound_id = 3;
switch Compound_id
    case 1
        ONCellFileName = 'ONRGCs_Connect4_04302401_movespot04302401.mat';
        ON_Simulation_name = 'trajectory04302401';
        compound_save_file_name ='ONconnect4_01.mat';
    case 2
        ONCellFileName = 'ONRGCs_Connect1_04302402_movespot04302402.mat';
        ON_Simulation_name = 'trajectory04302402';
        compound_save_file_name ='ONconnect1_01.mat';
    case 3
        time_stamp = '04302403';
        ONCellFileName = sprintf('ONRGCs_Connect2_%s_movespot%s.mat', time_stamp, time_stamp);
        ON_Simulation_name = sprintf('trajectory%s', time_stamp);
        compound_save_file_name ='ONconnect2_01.mat';
end
datafolder = './Simulation/Data';
is_trajectory = 1;

%% load OFF responses
loaddata = load(fullfile(datafolder, OFFCellFileName), 'num_RGCs', 'simulation_num_frame', 'RGC2NF', 'CellType',...
    'sim_id', 'pos_NFs', 'rf_dia', 'num_NFs', 'spot_rad');
num_RGCs = loaddata.num_RGCs;
simulation_num_frame = loaddata.simulation_num_frame;
RGC2NF_OFF = loaddata.RGC2NF;
CellType = loaddata.CellType;
sim_id = loaddata.sim_id;
pos_NFs = loaddata.pos_NFs;
rf_dia = loaddata.rf_dia;
num_NFs = loaddata.num_NFs;
spot_rad = loaddata.spot_rad;
OFF_Simulation_name = 'trajectory05';
assert(CellType == 0);
%%  Vet traces
sim_traces_OFF = nan(num_RGCs, simulation_num_frame);
foldername = sprintf('./Simulation/RFs/%s', sim_id);
for i = 1:num_RGCs
    rf_filename = sprintf('%s_%d.mat', OFF_Simulation_name, i);
    load(fullfile(foldername, rf_filename), 'simulation_trace');
    switch CellType
        case 0
            a = simulation_trace+0.5;
        case 1
            a = simulation_trace;
    end
    a(a<0) = 0;
    sim_traces_OFF(i, :) = a;
end
sim_traces_OFF = sim_traces_OFF-sim_traces_OFF(:, 1);

%% load ON responses
loaddata = load(fullfile(datafolder, ONCellFileName), 'num_RGCs', 'simulation_num_frame', 'CellType',...
    'sim_id', 'pos', 'pixel2um', 'square_length', 'extending_length', 'Fz', 'trajectoryfolder', 'trajectory_name',...
    'sim_fraction', 'num_connect');
num_RGCs = loaddata.num_RGCs;
simulation_num_frame = loaddata.simulation_num_frame;
CellType = loaddata.CellType;
sim_id = loaddata.sim_id;
pos = loaddata.pos;
square_length = loaddata.square_length;
extending_length = loaddata.extending_length;
pixel2um = loaddata.pixel2um;
Fz = loaddata.Fz;
trajectoryfolder = loaddata.trajectoryfolder;
trajectory_name = loaddata.trajectory_name;
sim_fraction = loaddata.sim_fraction;
num_connect = loaddata.num_connect;
assert(CellType == 1);

% create NF that corresponding to the OFF
RGC2NF_ON = distanceConnection(pos_NFs, pos, rf_dia);

%%  Vet traces
sim_traces_ON = nan(num_RGCs, simulation_num_frame);
foldername = sprintf('./Simulation/RFs/%s', sim_id);
for i = 1:num_RGCs
    rf_filename = sprintf('%s_%d.mat', ON_Simulation_name, i);
    load(fullfile(foldername, rf_filename), 'simulation_trace');
    switch CellType
        case 0
            a = simulation_trace+0.5;
        case 1
            a = simulation_trace;
    end
    a(a<0) = 0;
    sim_traces_ON(i, :) = a;
end
sim_traces_ON = sim_traces_ON-sim_traces_ON(:, 1);
%% Get NF responses
sim_NF = RGC2NF_OFF*sim_traces_OFF + RGC2NF_ON*sim_traces_ON;
% sim_NF = RGC2NF_OFF*sim_traces_OFF;

figure;
imagesc(sim_NF); colorbar

%%
spot_pos = load(fullfile(trajectoryfolder, trajectory_name), 'num_dots', 'positions');
num_dots = round(spot_pos.num_dots*sim_fraction);
spot_pos = spot_pos.positions;
spot_pos = spot_pos(1:num_dots, :);
%% Visualize traces
close all
is_generate_video = 1;
display_rad = 200;
display_rad = display_rad/pixel2um;
num_color = 256;
Colors = parula(num_color);
val_range = [min(sim_NF(:)) max(sim_NF(:))];
t = (0:simulation_num_frame-1)/Fz;
display_scalar = 200;
%
if is_generate_video
    moviename = sprintf('%s_%s_NF.mp4', sim_id, ON_Simulation_name);
    moviefolder = './Simulation/Movies';
    v = VideoWriter(fullfile(moviefolder, moviename),'MPEG-4');
    v.Quality = 95;
    v.FrameRate = 10;
    open(v)
end
%%

figure;
rectangle('Position',[-0.5*square_length-extending_length -0.5*square_length-extending_length,...
    square_length + extending_length*2 square_length + extending_length*2], 'EdgeColor', 'k');
keyboard;
%%
for i = 1:simulation_num_frame
    clf
    hold on
    for j = 1:num_NFs
        cval = sim_NF(j, i);
        cid = round(num_color*(cval-(val_range(1)))/range(val_range) + 1);
        if cid < 1
            cid = 1;
        elseif cid > num_color
            cid = num_color;
        end
        viscircles(pos_NFs(j, :), rf_dia/2, 'Color',Colors(cid, :));
    end

    if is_trajectory
        viscircles(spot_pos(i, :), spot_rad, 'Color',zeros(1, 3));
    else
        viscircles([(-0.5*square_length - extending_length) + spot_spd*i/Fz 0], spot_rad, 'Color',zeros(1, 3));
    end
    rectangle('Position',[-0.5*square_length-extending_length -0.5*square_length-extending_length,...
        square_length + extending_length*2 square_length + extending_length*2], 'EdgeColor', 'k');
    text(0.5*square_length+extending_length-330, 0.5*square_length+extending_length+50, sprintf('%.02G (s)', t(i)));
    if i == 1
        cmovdistance = spot_pos(i, :) - 0;
    else
        cmovdistance = spot_pos(i, :) - spot_pos(i-1, :);
    end
    if is_trajectory
        text(0.5*square_length+extending_length-170, 0.5*square_length+extending_length+50, sprintf('%d (um/s)',...
            round(sqrt(sum((cmovdistance*Fz).^2))*pixel2um)));
    else
        text(0.5*square_length+extending_length-170, 0.5*square_length+extending_length+50, sprintf('%d (um/s)',...
            round(spot_spd*pixel2um)));
    end
    xlim([-0.5*square_length-extending_length 0.5*square_length+extending_length]);
    ylim([-0.5*square_length-extending_length 0.5*square_length+extending_length]);
    plot(-0.5*square_length-extending_length+[100 100+display_scalar/pixel2um],...
        -(0.5*square_length+extending_length-50)*ones(1, 2), 'k', 'LineWidth', 2);
    plot((-0.5*square_length-extending_length+100)*ones(1, 2),...
        -(0.5*square_length+extending_length-50)+[0 0+display_scalar/pixel2um], 'k', 'LineWidth', 2);
    text(-0.5*square_length-extending_length+120, -0.5*square_length-extending_length+100,...
        sprintf('%d (um)', display_scalar));
    box off
    axis off
    pause(0.2)
    hold off
    if is_generate_video
        frame = getframe(gcf);
        writeVideo(v,frame);
    end
end
if is_generate_video
    close(v);
end

%%
compound_datafolder = './Simulation/Data/Compound';
save(fullfile(compound_datafolder, compound_save_file_name), 'sim_NF', 'pos_NFs', 'num_connect');
%%
[coeff,score,latent,tsquared,explained,mu] = pca(sim_NF');
%%
figure;
subplot(1, 2, 1); hold on
plot(score(:, 1), 'k')
plot(score(:, 2), 'r')
subplot(1, 2, 2); hold on
plot(score(:, 1), score(:, 2), 'k');
%%
num_history = 5;
summary_data_history = nan(num_history, 3);
for i = 1:num_history
    fea_ids = find(explained>5);
    data = [score(:, fea_ids) spot_pos];
    n = length(fea_ids);
    m = i-1;
    [beta_x, beta_y, Y_pred] = predictXY(data, n, m);
    corr_x = corr(diff(spot_pos((m+1):end, 1)), diff(Y_pred(:, 1)));
    corr_y = corr(diff(spot_pos((m+1):end, 2)), diff(Y_pred(:, 2)));
    % distance
    distance_pred = mean(sqrt(sum((spot_pos((m+1):end, :) - Y_pred).^2)))*pixel2um;
    summary_data_history(i, :) = [corr_x, corr_y, distance_pred];
end

m = 0;
num_latent = 10;
summary_data_numlatent = nan(num_latent, 3);
for i = 1:num_latent
    fea_ids = 1:i;
    data = [score(:, fea_ids) spot_pos];
    n = length(fea_ids);
    [beta_x, beta_y, Y_pred] = predictXY(data, n, m);
    corr_x = corr(diff(spot_pos((m+1):end, 1)), diff(Y_pred(:, 1)));
    corr_y = corr(diff(spot_pos((m+1):end, 2)), diff(Y_pred(:, 2)));
    % distance
    distance_pred = mean(sqrt(sum((spot_pos((m+1):end, :) - Y_pred).^2)))*pixel2um;
    summary_data_numlatent(i, :) = [corr_x, corr_y, distance_pred];
end
%%
figure; 
subplot(1, 3, 1);hold on
x = (1:num_history)-1;
plot(x, summary_data_history(:, 1), 'b');
plot(x, summary_data_history(:, 2), 'Color', 0.3*[0 0 1]);
ylim([0 1]);
xlabel('Number of history');
ylabel('Prediction correlation coefficient');
legend({'x', 'y'})
subplot(1, 3, 2);
plot(explained(1:num_latent)/100, 'b');
box off
xlabel('Principal components');
ylabel('Explained variance');
xlim([1 num_latent]);
ylim([0 1]);
subplot(1, 3, 3);hold on
plot(summary_data_numlatent(:, 1), 'b');
plot(summary_data_numlatent(:, 2), 'Color', 0.3*[0 0 1]);
xlim([1 num_latent]);
ylim([0 1]);
xlabel('Principal components');
ylabel('Prediction correlation coefficient');

%%
n = 10;
m = 0;
fea_ids = 1:n;
data = [score(:, fea_ids) spot_pos];
[beta_x, beta_y, Y_pred] = predictXY(data, n, m);
test_spot_pos = spot_pos((m+1):end, :);

figure; 
subplot(1, 2, 1); hold on 
Colors = parula(size(test_spot_pos, 1));
for i = 1:size(test_spot_pos, 1)-1
    plot(test_spot_pos([i i+1], 1), test_spot_pos([i i+1], 2), 'Color', Colors(i, :));
end
xlim([-600 600])
ylim([-600 600])
title('Spot trajectory');
subplot(1, 2, 2); hold on 
for i = 1:size(test_spot_pos, 1)-1
    plot(Y_pred([i i+1], 1), Y_pred([i i+1], 2), 'Color', Colors(i, :));
end
xlim([-600 600])
ylim([-600 600])
title('NF predicted trajectory');
%%
manually_input_save_name = 'ON_4conn_OFFRGCs';
compound_datafolder = './Simulation/Data/Compound';
save(fullfile(compound_datafolder, manually_input_save_name), ...
    'n', 'm', 'summary_data_history',...
    'summary_data_numlatent', 'test_spot_pos', 'Y_pred');

