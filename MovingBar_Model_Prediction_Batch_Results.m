clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';

recording_sets = {'e100724', 'f100724', 'a101224', 'b101224', 'c101224',   'd101224', 'e101224',...
             'b101424', 'c101424', 'd101424', 'e101424', 'a101624',   'b101624', 'd101624', 'e101624',...
             'b101924', 'c101924', 'd101924', 'e101924', 'b103124',   'e103124', 'a110424',...
             'c110424',                       'f110424', 'g110424',   'a110924', 'b110924', 'c110924',...
             'a111224'};

Fz = 100;
num_recording = length(recording_sets);
all_corr = nan(num_recording, 8);
all_SC = nan(num_recording, 2);
is_plot = 0;

%%
process_version = 'GaussianFitting_processed_082025_1.mat';
folder_name = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite';
processedFile = fullfile(folder_name, process_version);

if exist(processedFile, 'file')
    WN = load(processedFile,'data_sets', 'cell_type', 'location');
    fprintf('Loaded fitted parameters for LN model\n');
else
    error('Fitted parameters not found. Run WhiteNoise_ONOFFalpha_Comparison.m first.');
end
%% run SpotSizeAnalysis_simple.m to get the following values
ON_Nasal_CSR = -0.006;
ON_Temporal_CSR = -0.09;
OFF_Nasal_CSR = -0.056;
OFF_Temporal_CSR = -0.062;
CSRStrength = 100;
%%
count_csr_type = nan(num_recording, 1);
%%
for ii = 1:num_recording
    recording_name = recording_sets{ii};
    PredictionResults = nan(1, 7);
    BaselineCorr = nan(1, 1);
    clear PredTraces
    
    % Find corresponding WN data for this recording
    wn_idx = find(strcmp(WN.data_sets, recording_name));
    if isempty(wn_idx)
        error('Recording %s not found in WN data, skipping CSR constraint', recording_name);
    else
        current_cell_type = WN.cell_type{wn_idx};
        current_location = WN.location{wn_idx};
        
        % Determine CSR value based on cell type and location
        if strcmp(current_cell_type, 'ON') && strcmp(current_location, 'Nasal')
            current_CSR = ON_Nasal_CSR;
            count_csr_type(ii) = 1;
        elseif strcmp(current_cell_type, 'ON') && strcmp(current_location, 'Temporal')
            current_CSR = ON_Temporal_CSR;
            count_csr_type(ii) = 2;
        elseif strcmp(current_cell_type, 'OFF') && strcmp(current_location, 'Nasal')
            current_CSR = OFF_Nasal_CSR;
            count_csr_type(ii) = 3;
        elseif strcmp(current_cell_type, 'OFF') && strcmp(current_location, 'Temporal')
            current_CSR = OFF_Temporal_CSR;
            count_csr_type(ii) = 4;
        else
            warning('Unknown cell type (%s) or location (%s) for %s, using default CSR', ...
                    current_cell_type, current_location, recording_name);
            current_CSR = 0.09; % default value
        end
    end
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
    
    % Sanity check: bar_type should match WN.cell_type
    if ~isempty(wn_idx)
        if ~strcmp(bar_type, current_cell_type)
            error('SANITY CHECK FAILED: bar_type (%s) does not match WN.cell_type (%s) for recording %s', ...
                    bar_type, current_cell_type, recording_name);
        else
            fprintf('Sanity check passed for %s: %s %s cell, using CSR = %.3f\n', ...
                    recording_name, current_cell_type, current_location, current_CSR);
        end
    end
    response_name = recording_name;
    load_recording_name = recording_name;
    
    % Pass CSR value to the fitting script
    if ~isempty(wn_idx)
        CSR_value = current_CSR;
        fprintf('Using CSR = %.3f for %s (%s %s)\n', CSR_value, recording_name, current_cell_type, current_location);
    else
        CSR_value = 0.09; % default value
        fprintf('Using default CSR = %.3f for %s (WN data not found)\n', CSR_value, recording_name);
    end
    
    MovingBar_LinearNL_Fitting
end