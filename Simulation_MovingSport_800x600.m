clear; close all; clc
%% Build RGCs and their connection
density_RGCs = 200;
extending_length = 100;
square_length = 1000; % in um
num_RGCs = square_length^2*density_RGCs/1e6;
num_connect = 4;
pos = tilingNeuronGenerator(num_RGCs, square_length + extending_length*2);
cluster_ids = fixedNumConnection(pos, num_connect);


%%
Colors = lines(size(cluster_ids, 1));

close all
figure;
subplot(1, 2, 1); hold on
scatter(pos(:, 1), pos(:, 2), 5, 'k', 'filled');

%%
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
%% Build Receptive field of NF neurons
num_NFs = 50;
pixel2um = 2.5;
rf_dia = 8*32.5/pixel2um;
pos_NFs = tilingNeuronGenerator(num_NFs, square_length)+extending_length;


subplot(1, 2, 2); hold on
scatter(pos(:, 1), pos(:, 2), 5, 'k', 'filled');
scatter(pos_NFs(:, 1), pos_NFs(:, 2), 5, 'r', 'filled');
viscircles(pos_NFs, rf_dia/2, 'Color','m');

%% connect NF to RGCs by distance based weighted connectivity
RGC2NF = distanceConnection(pos_NFs, pos, rf_dia);
figure;
imagesc(RGC2NF); colorbar

%%
load('./Results/SummaryData_042324.mat', 'DataInfo', 'Data_TF', 'data_optimal_params_1', 'X','Y');
cids = DataInfo(:, 1)*2+DataInfo(:, 2);
sum_params = mean(data_optimal_params_1(cids == 0, :), 1);
SpatialKernel = gaussian2d(X, Y, sum_params);


CellType = 1;
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
spot_rad = 200; % in um
pixel2um = 2.5;
sim_id = '04242401';
surround_rad_fac = 1.5;

spot_rad = spot_rad/pixel2um;
width = size(SpatialKernel, 2);
height = size(SpatialKernel, 1);
spot_locs = [-spot_rad spot_rad];
[X, Y] = meshgrid(1:width, 1:height);
num_frame = length(TemporalKernel);
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
    foldername = sprintf('./Simulation/RFs/%s', sim_id);
    if ~exist(foldername, 'dir')
        mkdir(foldername)
    end
    save(fullfile(foldername, rf_filename), 'RF', 'SpatialKernel', 'TemporalKernel');
    clc
    fprintf('create single RF %d/%d', i, num_RGCs);
end

%% Generate RF based on the connectivity
weight_range = [0.8 1];

for i = 1:num_RGCs
    cid = find(cluster_ids == i);
    if isempty(cid)
        continue
    end
    [rid, cid] = ind2sub(size(cluster_ids), cid);
    weights = rand(1, num_connect-1)*range(weight_range)+0.8;
    weights = [1 weights];
    weights = weights/sum(weights);
    cids = [cluster_ids(rid, cid) cluster_ids(rid, ~ismember(1:4, cid))];


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
    % if sum(SpatialKernel(:)) ~= 0
    %     figure;
    %     imagesc(SpatialKernel); colorbar;
    %     keyboard;
    % end
    rf_filename = sprintf('connected_%d.mat', i);
    RF = new_RF;
    clear new_RF
    save(fullfile(foldername, rf_filename), 'RF', 'SpatialKernel', 'TemporalKernel');
    clc
    fprintf('create connected RF of RGC %d/%d', i, num_RGCs);
end


%% Simulate responses of a moving spots
spot_spd = 8000; % in um/s
spot_spd = spot_spd/pixel2um;
Fz = 100;
Simulation_name = 'movespot01';
Simulation_Time = 1; % in second
simulation_num_frame = round(Simulation_Time*Fz);

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
        spot_center = spot_locs(1) + spot_spd*i/Fz;
        pids = sqrt((X-spot_center).^2 + (Y-300).^2) < spot_rad;
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
display_rad = 100;
display_rad = display_rad/pixel2um;
num_color = 256;
Colors = parula(num_color);
val_range = [min(sim_traces(:)) max(sim_traces(:))];
t = (0:simulation_num_frame-1)/Fz;
display_scalar = 200;
%
moviename = sprintf('%s_%s.mp4', sim_id, Simulation_name);
moviefolder = './Simulation/Movies';
v = VideoWriter(fullfile(moviefolder, moviename),'MPEG-4');
v.Quality = 95;
v.FrameRate = 10;
open(v)
%
figure;
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

    viscircles([spot_locs(1) + spot_spd*i/Fz 300], spot_rad, 'Color',zeros(1, 3));
    rectangle('Position',[0 0 800 600], 'EdgeColor', 'k');
    text(10, -70, sprintf('%.02G (s)', t(i)));
    xlim([0 square_length + extending_length*2]);
    ylim([0 square_length + extending_length*2]);
    plot([100 100+display_scalar/pixel2um], 1000*ones(1, 2), 'k', 'LineWidth', 2);
    plot(100*ones(1, 2), [1000 1000+display_scalar/pixel2um], 'k', 'LineWidth', 2);
    text(150, 1100, sprintf('%d (um)', display_scalar));
    box off
    axis off
    pause(0.2)
    hold off
    frame = getframe(gcf);
    writeVideo(v,frame)
end
close(v)

%% Get NF responses
sim_NF = RGC2NF*sim_traces;
figure;
imagesc(sim_NF); colorbar

%% Visualize traces
close all
is_generate_video = 1;
display_rad = 100;
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
%
figure;
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

    viscircles([spot_locs(1) + spot_spd*i/Fz 300], spot_rad, 'Color',zeros(1, 3));
    rectangle('Position',[0 0 800 600], 'EdgeColor', 'k');
    text(10, -70, sprintf('%.02G (s)', t(i)));
    xlim([0 square_length + extending_length*2]);
    ylim([0 square_length + extending_length*2]);
    plot([100 100+display_scalar/pixel2um], 1000*ones(1, 2), 'k', 'LineWidth', 2);
    plot(100*ones(1, 2), [1000 1000+display_scalar/pixel2um], 'k', 'LineWidth', 2);
    text(150, 1100, sprintf('%d (um)', display_scalar));
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



