close all; clear; clc;
%%
ON1 = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar\a082924_moving_bar_processed.mat');
ON2 = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar\c082924_moving_bar_processed.mat');
OFF1 = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar\e100724_moving_bar_processed.mat');
OFF2 = load('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingBar\f100724_moving_bar_processed.mat');
%%
Fz = 100;
disp_direction = 180;
disp_contrast = 0;
disp_bar_witdth = [50, 200, 800];
disp_speeds = [500, 2000, 8000];
max_t = 459;
ct = (0:max_t-1)/Fz;
Colors = parula(4);
for i = 1:length(disp_contrast)
    figure; 
    for j = 1:length(disp_bar_witdth)
        for q = 1:length(disp_speeds)
            subplot(length(disp_bar_witdth), length(disp_speeds), (j-1)*length(disp_speeds)+q); hold on
            dir_id = find(ON1.dim1_moving_direction == disp_direction);
            ctr_id = find(ON1.dim2_contrast == disp_contrast);
            bw_id = find(ON1.dim3_bar_width == disp_bar_witdth(j));
            sp_id = find(ON1.dim4_speeds == disp_speeds(q));
            csig = squeeze(mean(ON1.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
            plot(ct, csig(1:max_t), 'Color', 0.7*[1 0 0]);

            dir_id = find(ON2.dim1_moving_direction == disp_direction);
            ctr_id = find(ON2.dim2_contrast == disp_contrast);
            bw_id = find(ON2.dim3_bar_width == disp_bar_witdth(j));
            sp_id = find(ON2.dim4_speeds == disp_speeds(q));
            csig = squeeze(mean(ON2.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
            plot(ct, csig(1:max_t), 'Color', 0.4*[1 1 0]);

            dir_id = find(OFF1.dim1_moving_direction == disp_direction);
            ctr_id = find(OFF1.dim2_contrast == disp_contrast);
            bw_id = find(OFF1.dim3_bar_width == disp_bar_witdth(j));
            sp_id = find(OFF1.dim4_speeds == disp_speeds(q));
            csig = squeeze(mean(OFF1.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
            plot(ct, csig(1:max_t), 'Color', 0.1*[0 0 0]);

            dir_id = find(OFF2.dim1_moving_direction == disp_direction);
            ctr_id = find(OFF2.dim2_contrast == disp_contrast);
            bw_id = find(OFF2.dim3_bar_width == disp_bar_witdth(j));
            sp_id = find(OFF2.dim4_speeds == disp_speeds(q));
            csig = squeeze(mean(OFF2.Data(dir_id, ctr_id, bw_id, sp_id, :, :), 5));
            plot(ct, csig(1:max_t), 'Color', 0.5*[0 0 1]);
            ylim([0 250]);
            xlim([ct(1) ct(end)]);
            xlabel('Time (s)');
            ylabel('Firing rate (spike/s)')
        end
    end
    sgtitle(sprintf('Direction: %d  Contrast: %0.2G', disp_direction, 1-disp_contrast));
end
