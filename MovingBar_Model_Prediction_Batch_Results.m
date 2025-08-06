clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_sets = {'e100724', 'f100724', 'a101224', 'b101224', 'c101224', 'd101224', 'e101224',...
             'b101424', 'c101424', 'd101424', 'e101424', 'a101624', 'b101624', 'd101624', 'e101624',...
             'b101924', 'c101924', 'd101924', 'e101924'};
% recording_sets = {'b101924'};
Fz = 100;
implement_case_id = 6;
num_recording = length(recording_sets);

for ii = 1:num_recording
    recording_name = recording_sets{ii};
    PredictionResults = nan(2, 2);
    BaselineCorr = nan(2, 1);
    clear PredTraces
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
    end
    response_name = recording_name;
    load_recording_name = recording_name;
    for jj = 1:2
        switch jj
            case 1
                is_blurry = 0;
            case 2
                is_blurry = 1;
        end
        MovingBar_LinearNL_Fitting
    end
    
end