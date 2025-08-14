clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
% recording_sets = {'b101424', 'c101424', 'd101424', 'e101424', 'a101624', 'b101624', 'd101624', 'e101624', 'd101924', 'e101924'};
recording_sets = {'b103124',   'e103124', 'a110424',...
             'c110424',                       'f110424', 'g110424',   'a110924', 'b110924', 'c110924',...
             'a111224'};
Fz = 100;
implement_case_id = 6;
num_recording = length(recording_sets);
is_only_NL = 0;
for ii = 1:num_recording
    recording_name = recording_sets{ii};
    
    switch recording_name
        case 'e100724'
            stim_wn_id = '003';
            stim_mb_id = '005';
            bar_type = 'OFF';
        case 'f100724'
            stim_wn_id = '001';
            stim_mb_id = '003';
            bar_type = 'OFF';
        case 'a101224'
            stim_wn_id = '002';
            stim_mb_id = '004';
            bar_type = 'OFF';
        case 'b101224'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'ON';
        case 'c101224'
            stim_wn_id = '001';
            stim_mb_id = '003';
            bar_type = 'OFF';
        case 'd101224'
            stim_wn_id = '001';
            stim_mb_id = '003';
            bar_type = 'ON';
        case 'e101224'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'OFF';
       case 'b101424'
            stim_wn_id = '002';
            stim_mb_id = '004';
            bar_type = 'OFF';
       case 'c101424'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'OFF';
       case 'd101424'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'ON';
       case 'e101424'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'ON';
       case 'a101624'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'ON';
       case 'b101624'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'ON';
       case 'd101624'
            stim_wn_id = '001';
            stim_mb_id = '003';
            bar_type = 'ON';
       case 'e101624'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'ON';
        case 'b101924'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'ON';
        case 'c101924'
            stim_wn_id = '002';
            stim_mb_id = '004';
            bar_type = 'OFF';
        case 'd101924'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'OFF';
        case 'e101924'
            stim_wn_id = '002';
            stim_mb_id = '004';
            bar_type = 'OFF';
         case 'b103124'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'ON';
        case 'e103124'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'OFF';
        case 'a110424'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'ON';
        case 'c110424'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'ON';
        case 'f110424'
            stim_wn_id = '002';
            stim_mb_id = '003';
            bar_type = 'OFF';
        case 'g110424'
            stim_wn_id = '001';
            stim_mb_id = '002';
            bar_type = 'OFF';
        case 'a110924'
            stim_wn_id = '004';
            stim_mb_id = '005';
            bar_type = 'ON';
        case 'b110924'
            stim_wn_id = '008';
            stim_mb_id = '009';
            bar_type = 'OFF';
        case 'c110924'
            stim_wn_id = '003';
            stim_mb_id = '005';
            bar_type = 'OFF';
        case 'a111224'
            stim_wn_id = '003';
            stim_mb_id = '004';
            bar_type = 'ON';
    end
    response_name = recording_name;
    load_recording_name = recording_name;
    Get_LinearNL_Params
    if ~ is_only_NL
        for jj = 1:2
            switch jj
                case 1
                    is_blurry = 0;
                case 2
                    is_blurry = 1;
            end
            MovingBar_LinearNL_Simulation
        end
    end
    
end
%%
keyboard;

%%
load(sprintf('./Results/MovingBar/%s', save_file_name), 'dim1_moving_direction', 'dim2_contrast',...
    'dim3_bar_width', 'dim4_speeds','dim5_repeats', 'dim6_time', 'resp', 'resp_s', 'cntr');

%%
dr_id = 3;
type_display = 4;
switch type_display
    case 0
        cresp = resp-0.3*resp_s;
    case 1
        cresp = resp;
    case 2
        cresp = resp-0.25*resp_s; 
    case 3
        cresp = resp-0.5*(resp_s-7e-5); 
    case 4
        cresp = resp-resp_s;
end
ct = (0:size(Data, 6)-1)/Fz;
close all
for q = 1:length(dim2_contrast)
    figure;
    for i = 1:length(dim3_bar_width)
        for j = 1:length(dim4_speeds)
            subplot(length(dim3_bar_width), length(dim4_speeds), (i-1)*length(dim4_speeds) + j);
            hold on
            dsig = squeeze(mean(Data(dr_id, q, i, j, :, :), 4, 'omitnan'));
            dsig = mean(dsig, 1);
            plot(ct, dsig, 'k');
            csig = squeeze(cresp(q, i, j, :));
            switch type_display
                case 0
                    csig = csig*max(dsig)./max(csig);
                    plot(ct, csig, 'm');

                case 1
                    plot(ct, 1.8*(nl_fuc(csig/divider)-15), 'm');
                    ylim([0 150]);
                case 2
                    plot(ct, 2.2*(nl_fuc(csig/divider)-25), 'm');
                    ylim([0 150]);
                case 3
                    plot(ct, 5*(nl_fuc(csig/divider)-55), 'm');
                    ylim([0 150]);
                case 4
                    plot(ct, 3.5*(nl_fuc(csig/divider)-35), 'm');
                    ylim([0 150]);
            end
            
            if i == length(dim3_bar_width)
                xlabel('Time (s)');
            elseif i == 1
                title(sprintf('%d (um / s)', dim4_speeds(j)));
            end
            if j == 1
                ylabel(sprintf('Bar width %d (um) \n Firing rate (spike/s)', dim3_bar_width(i)));
            end
            xlim([0 ct(end)]);
        end
    end
    sgtitle(sprintf('%s background contrast %.2G', recording_name, dim2_contrast(q)));

end

%%
exc = resp(1, 3, 1, :);
inh = resp_s(1, 3, 1, :);
cmb = resp(1, 3, 1, :)-(resp_s(1, 3, 1, :));
figure; hold on
plot(ct, squeeze(exc), 'k');
plot(ct, squeeze(inh), 'r');
plot(ct, squeeze(cmb), 'b');
title('Linear part');
xlabel('Time (s)')
legend({'Exc', 'Inh', 'Combined'})

