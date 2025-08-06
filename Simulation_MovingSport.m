clear; close all; clc
%% Task 
% implement speed bump
%% Build RGCs and their connection
time_stamp = '04302403';
CellType = 1;
num_connect = 2;

density_RGCs = 200;
extending_length = 100;
square_length = 1000; % in pixel
num_NFs = 70;
pixel2um = 2.5;
spot_rad = 200; % in um
surround_rad_fac = 1.5;
random_seeds = 100;

switch CellType
    case 1
        sim_id = sprintf('ONRGCs_Connect%d_%s', num_connect, time_stamp);
    case 0
        sim_id = sprintf('OFFRGCs_Connect%d_%s', num_connect, time_stamp);
end

rng(random_seeds)
num_RGCs = (square_length + extending_length*2)^2*density_RGCs/1e6;
pos = tilingNeuronGenerator(num_RGCs, square_length + extending_length*2);
pos = pos - 0.5*(square_length + extending_length*2);

cluster_ids = fixedNumConnection(pos, num_connect);


%%
Colors = lines(size(cluster_ids, 1));
close all
figure;
subplot(1, 2, 1); hold on
scatter(pos(:, 1), pos(:, 2), 5, 'k', 'filled');

%%
if num_connect > 1
    for i = 1:size(cluster_ids, 1)
        cids = cluster_ids(i, :);
        [~, mid] = min(sqrt(sum((pos(cids, :) - mean(pos(cids, :), 1)).^2, 2)));
        cpos = pos(cids, :);
        mpos = cpos(mid, :);
        cpos(mid, :) = [];
        for j = 1:(num_connect-1)
            plot([mpos(1) cpos(j, 1)], [mpos(2) cpos(j, 2)], 'Color', Colors(i, :));
        end
    end
end
%% Build Receptive field of NF neurons
rng(random_seeds)
rf_dia = 8*32.5/pixel2um;
pos_NFs = tilingNeuronGenerator(num_NFs, square_length);
pos_NFs = pos_NFs-0.5*square_length;

subplot(1, 2, 2); hold on
scatter(pos(:, 1), pos(:, 2), 5, 'k', 'filled');
scatter(pos_NFs(:, 1), pos_NFs(:, 2), 5, 'r', 'filled');
viscircles(pos_NFs, rf_dia/2, 'Color','m');

%% connect NF to RGCs by distance based weighted connectivity
RGC2NF = distanceConnection(pos_NFs, pos, rf_dia);
figure;
imagesc(RGC2NF); colorbar
keyboard;
%%
gird_axis = (-0.5*square_length - extending_length):(0.5*square_length + extending_length);
[X, Y] = meshgrid(gird_axis, gird_axis);
load('./Results/SummaryData_042324.mat', 'DataInfo', 'Data_TF', 'data_optimal_params_1');
cids = DataInfo(:, 1)*2+DataInfo(:, 2);
sum_params = mean(data_optimal_params_1(cids == 0, :), 1);
SpatialKernel = gaussian2d(X, Y, sum_params);


ca = Data_TF(DataInfo(:, 1) == CellType, :);
ca = ca./max(abs(ca), [], 2);

TemporalKernel = mean(ca, 1);
% SpatialKernel(SpatialKernel<0) = 0;
figure;
subplot(1, 2, 1)
imagesc(SpatialKernel); colorbar;
subplot(1, 2, 2)
plot(TemporalKernel)
%% Generate RF for each neurons
spot_rad = spot_rad/pixel2um;
% width = size(SpatialKernel, 2);
% height = size(SpatialKernel, 1);
spot_locs = [-spot_rad spot_rad];
% [X, Y] = meshgrid(1:width, 1:height);
num_frame = length(TemporalKernel);
foldername = sprintf('./Simulation/RFs/%s', sim_id);
if ~exist(foldername, 'dir')
    mkdir(foldername)
end
for i = 1:num_RGCs
    cparams = sum_params;
    cparams(1:2) = pos(i, :);
    SpatialKernel = gaussian2d(X, Y, cparams);
    mask = sqrt((X-pos(i, 1)).^2 + (Y-pos(i, 2)).^2) < spot_rad*surround_rad_fac;
    SpatialKernel(mask ~= 1) = 0;
    RF = repmat(SpatialKernel, 1, 1, num_frame).*reshape(TemporalKernel, 1, 1, []);

    switch CellType
        case 1
            RF = RF/(sum(RF(RF>0))+eps);
        case 0
            RF = RF/(abs(sum(RF(RF<0)))+eps);
    end
    if any(isnan(RF(:)))
        keyboard;
    end
    rf_filename = sprintf('single_%d.mat', i);
    
    save(fullfile(foldername, rf_filename), 'RF', 'SpatialKernel', 'TemporalKernel');
    clc
    fprintf('create single RF %d/%d', i, num_RGCs);
    close all
    figure; 
    imagesc(SpatialKernel); colorbar;
    pause(0.2);
end

%% Generate RF based on the connectivity
weight_range = [0.8 1];
is_display_RF = 0;
for i = 1:num_RGCs
    cid = find(cluster_ids == i);
    if isempty(cid)
        continue
    end
    if num_connect > 1
        [rid, cid] = ind2sub(size(cluster_ids), cid);
        weights = rand(1, num_connect-1)*range(weight_range)+0.8;
        weights = [1 weights];
        weights = weights/sum(weights);
        cids = [cluster_ids(rid, cid) cluster_ids(rid, ~ismember(1:num_connect, cid))];


        for j = 1:num_connect
            rf_filename = sprintf('single_%d.mat', cids(j));
            load(fullfile(foldername, rf_filename), 'RF', 'TemporalKernel');
            if j == 1
                new_RF = weights(j)*RF;
            else
                new_RF = new_RF + weights(j)*RF;
            end
        end
        SpatialKernel = squeeze(std(new_RF, [], 3));
        if is_display_RF
            if sum(SpatialKernel(:)) ~= 0
                figure;
                imagesc(SpatialKernel); colorbar;
                keyboard;
            end
        end
        rf_filename = sprintf('connected_%d.mat', i);
        RF = new_RF;
        clear new_RF
    else
        rf_filename = sprintf('single_%d.mat', cid);
        load(fullfile(foldername, rf_filename), 'RF', 'TemporalKernel');

    end
    save(fullfile(foldername, rf_filename), 'RF', 'SpatialKernel', 'TemporalKernel');
    clc
    fprintf('create connected RF of RGC %d/%d', i, num_RGCs);
end


%% Simulate responses of a moving spots
spot_spd = 8000; % in um/s
spot_spd = spot_spd/pixel2um;
Fz = 100;
is_trajectory = 0;
Simulation_name = sprintf('movespot%s', time_stamp);
Simulation_Time = 1; % in second
simulation_num_frame = round(Simulation_Time*Fz);
height = size(RF, 1);
width = size(RF, 2);
for j = 1:num_RGCs
    simulation_trace = nan(1,simulation_num_frame);
    cmov = ones(height, width, num_frame);
    rf_filename = sprintf('connected_%d.mat', j);
    if ~exist(fullfile(foldername, rf_filename), 'file')
        rf_filename = sprintf('single_%d.mat', j);
        load(fullfile(foldername, rf_filename), 'RF');
    else
        load(fullfile(foldername, rf_filename), 'RF');
    end
    for i = 1:simulation_num_frame
        simulation_trace(i) = cmov(:)'*RF(:);
        cmov(:, :, 1:end-1) = cmov(:, :, 2:end);
        spot_center = (-0.5*square_length - extending_length) + spot_spd*i/Fz;
        pids = sqrt((X-spot_center).^2 + (Y).^2) < spot_rad;
        blank = ones(height, width);
        blank(pids) = -1;
        cmov(:, :, end) = blank;
        % if i == 60
        %     figure;
        %     subplot(1, 3, 1);
        %     imagesc(SpatialKernel); colorbar;
        %     subplot(1, 3, 2);
        %     imagesc(blank, [-1 1]); colorbar;
        %     subplot(1, 3, 3);
        %     imagesc(squeeze(cmov(:, :, 1)), [-1 1])
        % end

        clc
        fprintf('simulation progress... %d/%d (%d/%d frame)', j, num_RGCs, i, simulation_num_frame);
    end
    rf_filename = sprintf('%s_%d.mat', Simulation_name, j);
    save(fullfile(foldername, rf_filename), 'simulation_trace');
end

%% Simulate responses of a moving spots in a trajactory
Simulation_name = sprintf('trajectory%s', time_stamp);
trajectory_name = 'trajectory_202404251111_s4.mat';
trajectoryfolder = './Simulation/Trajectories';
is_trajectory = 1;
sim_fraction = 0.25;
spot_pos = load(fullfile(trajectoryfolder, trajectory_name), 'num_dots', 'positions');
num_dots = round(spot_pos.num_dots*sim_fraction);
spot_pos = spot_pos.positions;
height = size(RF, 1);
width = size(RF, 2);
for j = 1:num_RGCs
    simulation_trace = nan(1,num_dots);
    cmov = ones(height, width, num_frame);
    rf_filename = sprintf('connected_%d.mat', j);
    if ~exist(fullfile(foldername, rf_filename), 'file')
        rf_filename = sprintf('single_%d.mat', j);
        load(fullfile(foldername, rf_filename), 'RF');
    else
        load(fullfile(foldername, rf_filename), 'RF');
    end
    for i = 1:num_dots
        simulation_trace(i) = cmov(:)'*RF(:);
        cmov(:, :, 1:end-1) = cmov(:, :, 2:end);
        pids = sqrt((X-spot_pos(i, 1)).^2 + (Y-spot_pos(i, 2)).^2) < spot_rad;
        blank = ones(height, width);
        blank(pids) = -1;
        cmov(:, :, end) = blank;
        % if i == 60
        %     figure;
        %     subplot(1, 3, 1);
        %     imagesc(SpatialKernel); colorbar;
        %     subplot(1, 3, 2);
        %     imagesc(blank, [-1 1]); colorbar;
        %     subplot(1, 3, 3);
        %     imagesc(squeeze(cmov(:, :, 1)), [-1 1])
        % end

        clc
        fprintf('simulation progress... %d/%d (%d/%d frame)', j, num_RGCs, i, simulation_num_frame);
    end
    rf_filename = sprintf('%s_%d.mat', Simulation_name, j);
    save(fullfile(foldername, rf_filename), 'simulation_trace');
end

simulation_num_frame =  num_dots;
%%
keyboard;
%%  Vet traces
sim_traces = nan(num_RGCs, simulation_num_frame);
for i = 1:num_RGCs
    rf_filename = sprintf('%s_%d.mat', Simulation_name, i);
    load(fullfile(foldername, rf_filename), 'simulation_trace');
    switch CellType
        case 0
            a = simulation_trace+0.5;
        case 1
            a = simulation_trace;
    end
    a(a<0) = 0;
    sim_traces(i, :) = a;
end
sim_traces = sim_traces-sim_traces(:, 1);

%% Visualize traces
close all
is_generate_video = 1;
display_rad = 100;
display_rad = display_rad/pixel2um;
num_color = 256;
Colors = parula(num_color);
val_range = [min(sim_traces(:)) max(sim_traces(:))];
t = (0:simulation_num_frame-1)/Fz;
display_scalar = 200;
%
if is_generate_video
    moviename = sprintf('%s_%s_ONAlpha.mp4', sim_id, Simulation_name);
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
%%

for i = 1:simulation_num_frame
    clf
    hold on
    for j = 1:num_RGCs
        cval = sim_traces(j, i);
        cid = round(num_color*(cval-(val_range(1)))/range(val_range) + 1);
        if cid < 1
            cid = 1;
        elseif cid > num_color
            cid = num_color;
        end
        viscircles(pos(j, :), display_rad, 'Color',Colors(cid, :));
    end

    % viscircles([spot_locs(1) + spot_spd*i/Fz 0], spot_rad, 'Color',zeros(1, 3));
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
        writeVideo(v,frame)
    end
end
if is_generate_video
    close(v)
end

%% Get NF responses
switch CellType
    case 0
        sim_NF = RGC2NF*sim_traces;
    case 1
        sim_NF = -RGC2NF*sim_traces;
end
figure;
imagesc(sim_NF); colorbar

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
    moviename = sprintf('%s_%s_NF.mp4', sim_id, Simulation_name);
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

    viscircles([(-0.5*square_length - extending_length) + spot_spd*i/Fz 0], display_rad, 'Color',zeros(1, 3));
    % viscircles([spot_locs(1) + spot_spd*i/Fz 0], display_rad, 'Color',zeros(1, 3));
    % viscircles(spot_pos(i, :), display_rad, 'Color',zeros(1, 3));
    rectangle('Position',[-0.5*square_length-extending_length -0.5*square_length-extending_length,...
        square_length + extending_length*2 square_length + extending_length*2], 'EdgeColor', 'k');
    text(0.5*square_length+extending_length-330, 0.5*square_length+extending_length+50, sprintf('%.02G (s)', t(i)));
    if i == 1
        cmovdistance = spot_pos(i, :) - 0;
    else
        cmovdistance = spot_pos(i, :) - spot_pos(i-1, :);
    end
    text(0.5*square_length+extending_length-170, 0.5*square_length+extending_length+50, sprintf('%d (um/s)',...
        round(spot_spd*pixel2um)));
    % text(0.5*square_length+extending_length-170, 0.5*square_length+extending_length+50, sprintf('%d (um/s)',...
    %     round(sqrt(sum((cmovdistance*Fz).^2))*pixel2um)));
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
datafolder = './Simulation/Data';
dataname = sprintf('%s_%s.mat', sim_id, Simulation_name);
save(fullfile(datafolder, dataname));
%%
% [coeff,score,latent,tsquared,explained,mu] = pca(sim_traces');
[coeff,score,latent,tsquared,explained,mu] = pca(sim_NF');
%%
figure;
subplot(1, 2, 1); hold on
plot(score(:, 1), 'k')
plot(score(:, 2), 'r')
subplot(1, 2, 2); hold on
plot(score(:, 1), score(:, 2), 'k')


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
figure; 
subplot(1, 2, 1); hold on 
Colors = parula(size(spot_pos, 1));
for i = 1:size(spot_pos, 1)-1
    plot(spot_pos([i i+1], 1), spot_pos([i i+1], 2), 'Color', Colors(i, :));
end
xlim([-600 600])
ylim([-600 600])
title('Spot trajectory');
subplot(1, 2, 2); hold on 
for i = 1:size(spot_pos, 1)-1
    plot(Y_pred([i i+1], 1), Y_pred([i i+1], 2), 'Color', Colors(i, :));
end
xlim([-600 600])
ylim([-600 600])
title('NF predicted trajectory');


