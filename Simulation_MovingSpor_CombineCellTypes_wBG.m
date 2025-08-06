clear; close all; clc
%%
Compound_id = 17;
is_generate_video = 0;
switch Compound_id
    case 1
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_05142403_diffmotion05142403.mat';
        OFF_Simulation_name = 'diffmotion05142403';
        ONCellFileName = 'ONRGCs_Connect4_sub0_05142404_diffmotion05142404.mat';
        ON_Simulation_name = 'diffmotion05142404';
        compound_save_file_name ='ONconnect4_01_BG60.mat';
    case 2
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg50_06072403_diffmotion06072403.mat';
        OFF_Simulation_name = 'diffmotion06072403';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg50_06072401_diffmotion06072401.mat';
        ON_Simulation_name = 'diffmotion06072401';
        compound_save_file_name ='ONconnect1_0607_BG50.mat';
    case 3
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg50_06072403_diffmotion06072403.mat';
        OFF_Simulation_name = 'diffmotion06072403';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg50_06072402_diffmotion06072402.mat';
        ON_Simulation_name = 'diffmotion06072402';
        compound_save_file_name ='ONconnect4_0607_BG50.mat';
    case 4
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg50_06112401_diffmotion06112401.mat';
        OFF_Simulation_name = 'diffmotion06112401';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg50_06112402_diffmotion06112402.mat';
        ON_Simulation_name = 'diffmotion06112402';
        compound_save_file_name ='ONconnect1_0611_BG50.mat';
    case 5
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg50_06112401_diffmotion06112401.mat';
        OFF_Simulation_name = 'diffmotion06112401';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg50_06112403_diffmotion06112403.mat';
        ON_Simulation_name = 'diffmotion06112403';
        compound_save_file_name ='ONconnect4_0611_BG50.mat';
    case 6
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg100_06132403_diffmotion06132403.mat';
        OFF_Simulation_name = 'diffmotion06132403';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg100_06132401_diffmotion06132401.mat';
        ON_Simulation_name = 'diffmotion06132401';
        compound_save_file_name ='ONconnect1_0613_BG100.mat';
    case 7
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg100_06132403_diffmotion06132403.mat';
        OFF_Simulation_name = 'diffmotion06132403';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg100_06132402_diffmotion06132402.mat';
        ON_Simulation_name = 'diffmotion06132402';
        compound_save_file_name ='ONconnect4_0613_BG100.mat';
    case 8
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg150_06152401_diffmotion06152401.mat';
        OFF_Simulation_name = 'diffmotion06152401';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg150_06152402_diffmotion06152402.mat';
        ON_Simulation_name = 'diffmotion06152402';
        compound_save_file_name ='ONconnect1_0615_BG150.mat';
    case 9
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg150_06152401_diffmotion06152401.mat';
        OFF_Simulation_name = 'diffmotion06152401';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg150_06152403_diffmotion06152403.mat';
        ON_Simulation_name = 'diffmotion06152403';
        compound_save_file_name ='ONconnect4_0615_BG150.mat';
    case 10
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg150_06182401_diffmotion06182401.mat';
        OFF_Simulation_name = 'diffmotion06182401';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg150_06182402_diffmotion06182402.mat';
        ON_Simulation_name = 'diffmotion06182402';
        compound_save_file_name ='ONconnect1_0618_BG150.mat';
    case 11
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg150_06182401_diffmotion06182401.mat';
        OFF_Simulation_name = 'diffmotion06182401';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg150_06182403_diffmotion06182403.mat';
        ON_Simulation_name = 'diffmotion06182403';
        compound_save_file_name ='ONconnect4_0618_BG150.mat';
    case 12
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg200_06192401_diffmotion06192401.mat';
        OFF_Simulation_name = 'diffmotion06192401';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg200_06192402_diffmotion06192402.mat';
        ON_Simulation_name = 'diffmotion06192402';
        compound_save_file_name ='ONconnect1_0619_BG200.mat';
    case 13
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg200_06192401_diffmotion06192401.mat';
        OFF_Simulation_name = 'diffmotion06192401';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg200_06192403_diffmotion06192403.mat';
        ON_Simulation_name = 'diffmotion06192403';
        compound_save_file_name ='ONconnect4_0619_BG200.mat';
    case 14
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg100_06202401_diffmotion06202401.mat';
        OFF_Simulation_name = 'diffmotion06202401';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg100_06202402_diffmotion06202402.mat';
        ON_Simulation_name = 'diffmotion06202402';
        compound_save_file_name ='ONconnect1_0620_BG100.mat';
    case 15
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg100_06202401_diffmotion06202401.mat';
        OFF_Simulation_name = 'diffmotion06202401';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg100_06202403_diffmotion06202403.mat';
        ON_Simulation_name = 'diffmotion06202403';
        compound_save_file_name ='ONconnect4_0620_BG100.mat';
    case 16
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg200_06212401_diffmotion06212401.mat';
        OFF_Simulation_name = 'diffmotion06212401';
        ONCellFileName = 'ONRGCs_Connect1_sub1_bg200_06212402_diffmotion06212402.mat';
        ON_Simulation_name = 'diffmotion06212402';
        compound_save_file_name ='ONconnect1_0621_BG200.mat';
    case 17
        OFFCellFileName = 'OFFRGCs_Connect1_sub0_bg200_06212401_diffmotion06212401.mat';
        OFF_Simulation_name = 'diffmotion06212401';
        ONCellFileName = 'ONRGCs_Connect4_sub1_bg200_06212403_diffmotion06212403.mat';
        ON_Simulation_name = 'diffmotion06212403';
        compound_save_file_name ='ONconnect4_0621_BG200.mat';
end
datafolder = './Simulation/Data';
is_trajectory = 1;

%% load OFF responses
loaddata = load(fullfile(datafolder, OFFCellFileName), 'num_RGCs', 'simulation_num_frame', 'RGC2NF', 'CellType',...
    'sim_id', 'pos_NFs', 'rf_dia', 'num_NFs', 'spot_rad', 'pos', 'pixel2um');
num_RGCs = loaddata.num_RGCs;
simulation_num_frame = loaddata.simulation_num_frame;
RGC2NF_OFF = loaddata.RGC2NF;
CellType = loaddata.CellType;
sim_id = loaddata.sim_id;
pos_NFs = loaddata.pos_NFs;
rf_dia = loaddata.rf_dia;
num_NFs = loaddata.num_NFs;
spot_rad = loaddata.spot_rad;
pos_OFF = loaddata.pos;
assert(CellType == 0);
%%  Vet traces
sim_traces_OFF = nan(num_RGCs, simulation_num_frame);
foldername = sprintf('./Simulation/RFs/%s', sim_id);
for i = 1:num_RGCs
    rf_filename = sprintf('%s_%d.mat', OFF_Simulation_name, i);
    load(fullfile(foldername, rf_filename), 'simulation_trace');
    switch CellType
        case 0
            a = simulation_trace;
        case 1
            a = simulation_trace+0.3;
    end
    a(a<0) = 0;
    sim_traces_OFF(i, :) = a;
end
% sim_traces_OFF = sim_traces_OFF-sim_traces_OFF(:, 1);

%% load ON responses
loaddata = load(fullfile(datafolder, ONCellFileName), 'num_RGCs', 'simulation_num_frame', 'CellType',...
    'sim_id', 'pos', 'pixel2um', 'patch_length_w', 'patch_length_h', 'extending_length', 'Fz', 'trajectoryfolder', 'trajectory_name',...
    'canvas', 'num_connect', 'Stimulus_speed_bg', 'Stimulus_track_t', 'Stimulus_speed_ob', 'pos');
num_RGCs = loaddata.num_RGCs;
simulation_num_frame = loaddata.simulation_num_frame;
CellType = loaddata.CellType;
sim_id = loaddata.sim_id;
pos = loaddata.pos;
patch_length_w = loaddata.patch_length_w;
patch_length_h = loaddata.patch_length_h;
Stimulus_speed_bg = loaddata.Stimulus_speed_bg;
Stimulus_speed_ob = loaddata.Stimulus_speed_ob;
Stimulus_track_t = loaddata.Stimulus_track_t;
canvas = loaddata.canvas;
extending_length = loaddata.extending_length;
pixel2um = loaddata.pixel2um;
Fz = loaddata.Fz;
num_connect = loaddata.num_connect;
pos_ON = loaddata.pos;
assert(CellType == 1);
ch = size(canvas, 1);
cw = size(canvas, 2);


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
            a = simulation_trace;
        case 1
            a = simulation_trace+0.3;
    end
    a(a<0) = 0;
    sim_traces_ON(i, :) = a;
end
% sim_traces_ON = sim_traces_ON-sim_traces_ON(:, 1);
%% Get NF responses
sim_NF = RGC2NF_OFF*sim_traces_OFF + RGC2NF_ON*sim_traces_ON;
% sim_NF = RGC2NF_OFF*sim_traces_OFF;

figure;
imagesc(sim_NF); colorbar

%%
%% Visualize traces
close all

display_rad = 200;
display_rad = display_rad/pixel2um;
num_color = 256;
Colors = parula(num_color);
val_range = [min(sim_NF(:)) max(sim_NF(:))];
t = (0:simulation_num_frame-1)/Fz;
display_scalar = 200;
%
if is_generate_video
    moviename = sprintf('%s_%s_NF_BG60.mp4', sim_id, ON_Simulation_name);
    moviefolder = './Simulation/Movies';
    v = VideoWriter(fullfile(moviefolder, moviename),'MPEG-4');
    v.Quality = 95;
    v.FrameRate = 10;
    open(v)
end
%%
figure;
rectangle('Position',[-0.5*patch_length_w -0.5*patch_length_h, patch_length_w patch_length_h], 'EdgeColor', 'k');
if is_generate_video
    keyboard;
end
%%
sw = patch_length_w;
sh = patch_length_h;
spot_center_bg = [0 0];
spot_center_ob = [0 0];
correct_coordi = [400, 200];
spot_pos = nan(simulation_num_frame, 2);
bg_pos = nan(simulation_num_frame, 2);
t = 0;
for i = 1:simulation_num_frame
    clf
    hold on
    t = t + 1/Fz;
    pid = sum(Stimulus_track_t<t)+1;
    spot_center_bg = [0, spot_center_bg(2) + Stimulus_speed_bg(pid)/Fz];
    bg_pos(i, :) = spot_center_bg;
    spot_center_ob_pre = spot_center_ob;
    spot_center_ob = [spot_center_ob(1) + Stimulus_speed_ob(pid)/Fz, 0];
    spot_pos(i, :) = spot_center_ob;
    
    %pids = sqrt((X-spot_center_ob(2)).^2 + (Y-spot_center_ob(1)).^2) < spot_rad;
    x_range = round(ch*0.5-0.5*sh-spot_center_bg(1))+1;
    x_range = x_range:(x_range+sh-1);
    y_range = round(cw*0.5-0.5*sw-spot_center_bg(2))+1;
    y_range = y_range:(y_range+sw-1);
    newpage = canvas(x_range, y_range);
    % newpage(pids(extending_length:extending_length+sh-1, extending_length:extending_length+sw-1)) = -1;
    imagesc(newpage);colormap(gray);
    for j = 1:num_NFs
        cval = sim_NF(j, i);
        cid = round(num_color*(cval-(val_range(1)))/range(val_range) + 1);
        if cid < 1
            cid = 1;
        elseif cid > num_color
            cid = num_color;
        end
        viscircles(pos_NFs(j, :)+correct_coordi, rf_dia/2, 'Color',Colors(cid, :));
    end

    viscircles(spot_center_ob + correct_coordi, spot_rad, 'Color',zeros(1, 3));
    text(0.5*patch_length_w-130, 0.5*patch_length_h+220, sprintf('%.02G (s)', t));
    cmovdistance = spot_center_ob - spot_center_ob_pre;

    text(0.5*patch_length_w-370, 0.5*patch_length_h+220, sprintf('%d (um/s)',...
        round(sqrt(sum((cmovdistance*Fz).^2))*pixel2um)));

    xlim([0 patch_length_w]);
    ylim([0 patch_length_h]);
    % plot(-0.5*square_length-extending_length+[100 100+display_scalar/pixel2um],...
    %     -(0.5*square_length+extending_length-50)*ones(1, 2), 'k', 'LineWidth', 2);
    % plot((-0.5*square_length-extending_length+100)*ones(1, 2),...
    %     -(0.5*square_length+extending_length-50)+[0 0+display_scalar/pixel2um], 'k', 'LineWidth', 2);
    % text(-0.5*square_length-extending_length+120, -0.5*square_length-extending_length+100,...
    %     sprintf('%d (um)', display_scalar));
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
save(fullfile(compound_datafolder, compound_save_file_name), 'sim_NF', 'pos_NFs', 'num_connect', 'spot_pos', 'bg_pos', 'canvas',...
     'RGC2NF_OFF', 'sim_traces_OFF', 'RGC2NF_ON', 'sim_traces_ON', 'pos_OFF', 'pos_ON');
%%
keyboard;
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

%% Test the x coordinate prediction (x axis at time)
n = 10;
nhist = 3;
fea_ids = 1:n;
data = [score(:, fea_ids) spot_pos];


test_spot_pos = spot_pos(1:end, :);
figure;  hold on 
t = (0:size(test_spot_pos, 1)-1)/Fz;
plot(t, test_spot_pos(:, 1), 'Color', 'k');
Colors = [linspace(0.5, 1, nhist)' zeros(nhist, 2)];
for i = 1:nhist

    m = i-1;
    [beta_x, beta_y, Y_pred] = predictXY(data, n, m);
    plot(t(m+1:end), Y_pred(:, 1), 'Color', Colors(i, :));
end
xlabel('Time (s)');
ylabel('x position');
% legend({'Spot trajectory', 'NF predicted trajectory'});

%% Test the x coordinate prediction (x axis at time)
n = 10;
nhist = 3;
fea_ids = 1:n;
data = [score(:, fea_ids) bg_pos(:, 2)];


test_spot_pos = bg_pos(1:end, :);
figure;  hold on 
t = (0:size(test_spot_pos, 1)-1)/Fz;
plot(t, test_spot_pos(:, 2), 'Color', 'k');
Colors = [linspace(0.5, 1, nhist)' zeros(nhist, 2)];
for i = 1:nhist

    m = i-1;
    [beta_x, beta_y, Y_pred] = predictXY(data, n, m);
    plot(t(m+1:end), Y_pred(:, 2), 'Color', Colors(i, :));
end
xlabel('Time (s)');
ylabel('x position');
% legend({'Spot trajectory', 'NF predicted trajectory'});

%%
n = 10;
m = 0;
fea_ids = 1:n;
data = [score(:, fea_ids) spot_pos(:, 1) bg_pos(:, 2)];
[beta_x, beta_y, Y_pred] = predictXY(data, n, m);
test_spot_pos = [spot_pos((m+1):end, 1) bg_pos((m+1):end, 2)];

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
manually_input_save_name = sprintf('%s_trajactory_result', compound_save_file_name);
compound_datafolder = './Simulation/Data/Compound';
save(fullfile(compound_datafolder, manually_input_save_name), ...
    'n', 'm', 'summary_data_history',...
    'summary_data_numlatent', 'test_spot_pos', 'Y_pred');

