clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_name = 'c081224';
load_recording_name = [recording_name '01'];
Fz = 100;
%% get the linear and nonlinear parts
load([load_data_folder load_recording_name '.mat'], 'masked_STAmat', 'PBs', 'FRs', 'stdSTA');
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_MovingNoise_1.mat'], 'OLED');
pix2um = OLED.pixelSize;
smtstdSTA = medfilt2(stdSTA);
[x, y, dist2d] = peak_distance(smtstdSTA);
minD = 100;
ythr = quantile(y(x>minD), 0.99);
binary_image = smtstdSTA>ythr;
[largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
fprintf('%0.3G', sqrt(sum(largest_segment_mask(:))*4/pi)*pix2um);
%%
num_gauss = 1;
image = stdSTA';
initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];
initial_params = repmat(initial_params, num_gauss, 1);
initial_params(:, 1:2) = initial_params(:, 1:2) + 10*rand(num_gauss, 2);
initial_params = initial_params';
objective_function = @(params) 1-corr(image(:), reshape(gaussian_multi(params, image, num_gauss), [], 1));
options.MaxFunEvals = 600*length(initial_params(:));
[optimal_params,fval,exitflag,output]= fminsearch(objective_function, initial_params);

gaussian_model = gaussian_multi(optimal_params, image, num_gauss);
figure; 
subplot(1, 2, 1)
imagesc(smtstdSTA'); hold on
plot([400 400], [0 600], '--k');
plot([0 800], [300 300], '--k');
subplot(1, 2, 2)
imagesc(gaussian_model); hold on
plot([400 400], [0 600], '--k');
plot([0 800], [300 300], '--k');
% Spatial temporal filter is in "masked_STAmat"

% Nonlinear parts can be derived from "PBs", "FRs"
%%
keyboard;
%%
[nl_fuc, divider]  = getNonlinearFunc(PBs, FRs);
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
imagesc(masked');
subplot(1, 2, 2)
imagesc(masked_STAmat(:, :, 45)');

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
plot(t, nl_fuc(resp*0.1/divider), 'k');

%% Simulation (1) varied size spots
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
load([stim_data_folder 'Temporal_AlphaRGC_' recording_name '_001_Retina_1_MovingBar_1.mat'], 'MB_IN');

loadFileName = sprintf('%s_moving_bar_processed.mat', recording_name);
load(sprintf('./Results/MovingBar/%s', loadFileName), 'dim1_moving_direction', 'dim2_bar_width', 'dim3_speeds',...
    'dim4_repeats', 'dim5_time', 'Data', 'Moving_end_time_ids');

Fz = 100;
%%
dr_id = 3;
bw_id = 5;
sp_id = 1;
figure;
% ct = (0:round(mean(Moving_end_time_ids(dr_id, bw_id, sp_id, :), 4, 'omitnan'))-1)/Fz
ct = (0:size(Data, 5)-1)/Fz;
plot(ct, squeeze(mean(Data(dr_id, bw_id, sp_id, :, :), 4, 'omitnan')));
title(sprintf('%d (um) %d (um/s)', dim2_bar_width(bw_id), dim3_speeds(sp_id)));

close all
cmov = -1*ones(size(masked_STAmat));
num_time = length(x);
t = (0:num_time-1)/Fz;
resp = nan(length(dim2_bar_width), length(dim3_speeds), size(Data, 5));
dim1 = size(cmov, 1);
dim2 = size(cmov, 2);

for k = 1:length(dim2_bar_width)
    for j = 1:length(dim3_speeds)
        bw = dim2_bar_width(k)/pix2um;
        cAng = dim1_moving_direction(dr_id);
        step = (dim3_speeds(j)/pix2um).*[cos(a2d(cAng+180)) sin(a2d(cAng+180))];
        posi = [cos(a2d(cAng))*0.5*diag sin(a2d(cAng))*0.5*diag]+...
            [cos(a2d(cAng))*0.5*bw sin(a2d(cAng))*0.5*bw];
        num_time = round(diag*Fz/step(1));
        csig = nan(1, num_time);
        for i = 1:num_time
            dstep = i*step/Fz;
            move = posi+dstep;
            switch cAng
                case 0
                    masked = abs((Y-0.5*(size(masked_STAmat, 1)+1)-move(1))) < 0.5*bw;
                case 180
                    masked = abs((Y-0.5*(size(masked_STAmat, 1)+1)-move(1))) < 0.5*bw;
            end
            csig(i) = mean(cmov.*masked_STAmat, 'all');
            cmov(:, :, 1) = [];
            canvas = 2*(double(masked)-0.5);
            cmov(:, :, end+1) = canvas;
            if mod(i, 10) == 1
                figure(1);
                imagesc(canvas', [-1 1]); colorbar;
            end
            clc
            fprintf('Progress ... %d/%d, %d/%d, %d/%d \n',  k, length(dim2_bar_width), j, length(dim3_speeds), i, num_time);
        end
        resp(k, j, 1:length(csig)) = csig;
    end
end

%% save
save_file_name = sprintf('%s_moving_bar_simulated.mat', recording_name);
save(sprintf('./Results/MovingBar/%s', save_file_name), 'dim1_moving_direction', 'dim2_bar_width', 'dim3_speeds',...
    'dim4_repeats', 'dim5_time', 'resp');

%%
ct = (0:size(Data, 5)-1)/Fz;
figure;
for i = 1:length(dim2_bar_width)
    for j = 1:length(dim3_speeds)
        subplot(length(dim2_bar_width), length(dim3_speeds), (i-1)*length(dim3_speeds) + j);
        hold on
        plot(ct, squeeze(mean(Data(dr_id, i, j, :, :), 4, 'omitnan')));
        plot(ct, nl_fuc(squeeze(resp(i, j, :))/divider));
        if i == length(dim2_bar_width)
            xlabel('Time (s)');
        elseif i == 1
            title(sprintf('%d (um / s)', dim3_speeds(j)));
        end
        ylim([0 150]);
        if j == 1
            ylabel(sprintf('Bar width %d (um) \n Firing rate (spike/s)', dim2_bar_width(i)));
        end
        xlim([0 ct(end)]);
    end
end
sgtitle(sprintf('%s', recording_name));

