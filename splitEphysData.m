function splitEphysData(matFilePath, ephysFilePath, minSectionDuration, samplingRate)

    if nargin < 3 || isempty(minSectionDuration)
        minSectionDuration = 1 * 60; % Default minimum section duration in seconds
    end
    if nargin < 4 || isempty(samplingRate)
        samplingRate = 10000; % Default sampling rate
    end

    % Load the ephys data structure
    if exist(matFilePath, 'file')
        load(matFilePath, 'ephysData');
    else
        error('File does not exist.');
    end

    detect_threshold = 2.5;
    % Process each entry in the ephysData struct
    FolderPath = [ephysFilePath 'sections\'];
    StimulusPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation';
    for i = 1:length(ephysData)
        % Check if the data has already been split
        dataName = ephysData(i).ephys_data_name;
        FileNames = dir(fullfile(FolderPath, '*.mat'));
        FileIds = find(contains({FileNames.name}, dataName));
        nFile = length(FileIds);
        alreadySplit = true;
        if nFile <= 0
            alreadySplit = false;
        end
        
        if alreadySplit
            disp(['Data for ' dataName ' has already been split.']);
            continue;
        end
        
        % Load the specific data file
        if exist([ephysFilePath dataName '.mat'], 'file')
            load([ephysFilePath dataName '.mat']);
        else
            disp(['Data file for ' dataName ' does not exist.']);
            continue;
        end
        
        % Identify the channel for processing
        selectedChannel = ephysData(i).ephys_channel;
        
        % Assuming data is in a variable named `data` within the .mat file
        triggerSignal = data(:, end);
        ephysSignal = data(:, selectedChannel);
        
         % Initialize variables for section detection
        sectionStarts = [];
        sectionEnds = [];
        minSamples = minSectionDuration * samplingRate; % Convert duration to samples
        
        
        % [problem sometimes the off part is much slower
        trigger_edges = tiggerdetection_continous(triggerSignal, detect_threshold);
        %%
%         trigger_edges_start = find(diff([triggerSignal(1); triggerSignal(:)]) > 2.4);
%         trigger_edges_end = find(diff([triggerSignal(1); triggerSignal(:)]) < -2.4);
%         
%         if trigger_edges_end(1) < trigger_edges_start(1)
%             trigger_edges_end(1) = [];
%         end
%         assert(trigger_edges_end(1) > trigger_edges_start(1))
%         trigger_edges = sort([trigger_edges_start; trigger_edges_end]);
%         if mod(length(trigger_edges), 2)==1
%             trigger_edges(end) = [];
%         end
%         trigger_edges = reshape(trigger_edges, 2, [])';
        trigger_edges(diff(trigger_edges, [], 2)<minSamples, :) = [];
        FileNames = dir(fullfile(StimulusPath, '*.mat'));
        FileIds = find(contains({FileNames.name}, ephysData(i).stimulation_name));
        numSections = length(FileIds);
        % Save each section into a separate file
        for sec = 1:numSections
            sectionData = ephysSignal(trigger_edges(sec, 1):trigger_edges(sec, 2));
            newFileName = [dataName, '_', num2str(sec), '.mat'];
            save([ephysFilePath 'sections/' newFileName], 'sectionData');
            fprintf('save ... %s \n', newFileName);
        end
    end
    disp('Data sections saved successfully.');
end
