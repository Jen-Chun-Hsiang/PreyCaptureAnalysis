clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
% recording_sets = {'b101424', 'c101424', 'd101424', 'e101424', 'a101624', 'b101624', 'd101624', 'e101624', 'd101924', 'e101924'};
% recording_sets = {'b103124',   'e103124', 'a110424',...
%              'c110424',                       'f110424', 'g110424',   'a110924', 'b110924', 'c110924',...
%              'a111224'};

recording_sets = {'e100724', 'f100724', 'a101224', 'b101224', 'c101224',   'd101224', 'e101224',...
             'b101424', 'c101424', 'd101424', 'e101424', 'a101624',   'b101624', 'd101624', 'e101624',...
             'b101924', 'c101924', 'd101924', 'e101924', 'b103124',   'e103124', 'a110424',...
             'c110424',                       'f110424', 'g110424',   'a110924', 'b110924', 'c110924',...
             'a111224'};
Fz = 100;
implement_case_id = 6;
num_recording = length(recording_sets);
is_only_NL = 0;
is_display = 0;
surround_sf_type = 'fixed'; % 'fixed', 'linearfit'

process_version = 'GaussianFitting_processed_082025_1.mat';
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
processedFile = fullfile(folder_name, process_version);

if exist(processedFile, 'file')
    load(processedFile, 'gauss_est', 'Gauss_TF_est', 'data_sets');
    fprintf('Loaded fitted parameters for LN model\n');
else
    error('Fitted parameters not found. Run WhiteNoise_ONOFFalpha_Comparison.m first.');
end
assert(length(data_sets) == size(gauss_est, 1), 'Data sets and fitted parameters do not match.');
assert(length(data_sets) == size(Gauss_TF_est, 1), 'Data sets and fitted parameters do not match.');

recording_name = recording_sets{1};
load([load_data_folder recording_name '.mat'], 'masked_STAmat', 'stdSTA');
[D1_mat, D2_mat, D3_mat] = size(masked_STAmat);
assert(size(stdSTA, 1) == D1_mat && size(stdSTA, 2) == D2_mat, 'stdSTA size mismatch');
clear masked_STAmat;

for ii = 1:num_recording   %num_recording
    recording_name = recording_sets{ii};
    

    loadFileName = sprintf('%s_moving_bar_processed.mat', recording_name);
    load(sprintf('./Results/MovingBar/%s', loadFileName), 'dim1_moving_direction', 'dim2_contrast', 'dim3_bar_width',...
    'dim4_speeds', 'dim5_repeats', 'dim6_time', 'Data', 'Moving_end_time_ids');
    cell_idx = find(strcmp(recording_name, data_sets));
    fprintf('Processing %s (%d/%d)...\n', recording_name, ii, num_recording);
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
            bar_type = 'ON';
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
    MovingBar_LinearNL_Simulation_steamline
end


