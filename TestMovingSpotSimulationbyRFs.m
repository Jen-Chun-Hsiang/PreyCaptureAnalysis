close all

%%
% load data from @DataSummaryPlot
%%
SpatialKernel = Z;
CellType = 1;

%%
switch CellType
    case 0
        tit_name = 'OFF cell';
    case 1
        tit_name = 'ON cell';
end
ca = Data_TF(DataInfo(:, 1) == CellType, :);
ca = ca./max(abs(ca), [], 2);


TemporalKernel = mean(ca, 1);
% SpatialKernel(SpatialKernel<0) = 0;
figure;
t = WinT(1):1/Fz:WinT(end);
plot(t(2:end), TemporalKernel)
xlabel('Times (s)')
box off
title(tit_name);


spot_rad = 200; % in um
spot_spd = 8000; % in um/s
pixel2um = 2.5;

spot_rad = spot_rad/pixel2um;
spot_spd = spot_spd/pixel2um;

width = size(SpatialKernel, 2);
height = size(SpatialKernel, 1);
spot_locs = [-spot_rad spot_rad];
[X, Y] = meshgrid(1:width, 1:height);
num_frame = length(TemporalKernel);
mask = sqrt((X-400).^2 + (Y-300).^2) < spot_rad*1.5;
SpatialKernel(mask ~= 1) = 0;
RF = repmat(SpatialKernel, 1, 1, num_frame).*reshape(TemporalKernel, 1, 1, []);
switch CellType
    case 1
        RF = RF/sum(RF(RF>0));
    case 0
        RF = RF/abs(sum(RF(RF<0)));
end

%%
Simulation_Time = 1; % in second
cmov = ones(height, width, num_frame);
simulation_num_frame = round(Simulation_Time*Fz);
simulation_trace = nan(1,simulation_num_frame);
for i = 1:simulation_num_frame
    rec_cmov = cmov(:);
    switch CellType
        case 0
            rec_cmov(rec_cmov>0) = 0;
        case 1
            rec_cmov(rec_cmov<0) = 0;
    end
    simulation_trace(i) = rec_cmov'*RF(:);
    cmov(:, :, 1:end-1) = cmov(:, :, 2:end);
    spot_center = spot_locs(1) + spot_spd*i/Fz;
    pids = sqrt((X-spot_center).^2 + (Y-300).^2) < spot_rad;
    blank = ones(height, width);
    blank(pids) = -1;
    cmov(:, :, end) = blank;
    if i == 60
        figure;
        subplot(1, 3, 1);
        imagesc(SpatialKernel); colorbar;
        subplot(1, 3, 2);
        imagesc(blank, [-1 1]); colorbar;
        subplot(1, 3, 3);
        imagesc(squeeze(cmov(:, :, 1)), [-1 1])
    end
    clc
    fprintf('simulation progress... %d/%d', i, simulation_num_frame);
end

%%
close all
switch CellType
    case 0
        a = simulation_trace;
    case 1
        a = simulation_trace-0.3;
end
t =(0:length(a)-1)/Fz;
figure; hold on
plot(t, a, 'k');
a(a<0) = 0;
plot(t, a, 'b');
xlabel('Time (s)');
ylabel('Effective contrast');
title(tit_name);


