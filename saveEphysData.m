function saveEphysData(matFilePath)
    % Check if the MAT file exists
    if exist(matFilePath, 'file')
        % Load the data from the MAT file
        load(matFilePath, 'ephysData');
    else
        % Initialize ephysData if the file doesn't exist
        ephysData = struct('ephys_data_name', {}, 'ephys_channel', {}, 'stimulation_name', {});
    end

    % Get user input
    ephys_data_name = input('Enter ephys data name (e.g., b120723_002): ', 's');
    ephys_channel = input('Enter ephys channel number: ');
    stimulation_name = input('Enter stimulation name (e.g., Temporal_AlphaRGC_c030124_001): ', 's');

    % Check if the input combination is unique
    isUnique = true;
    for i = 1:length(ephysData)
        if strcmp(ephysData(i).ephys_data_name, ephys_data_name) && ...
           ephysData(i).ephys_channel == ephys_channel && ...
           strcmp(ephysData(i).stimulation_name, stimulation_name)
            isUnique = false;
            break;
        end
    end

    % Add the new data if it's unique
    if isUnique
        newEntry = struct('ephys_data_name', ephys_data_name, ...
                          'ephys_channel', ephys_channel, ...
                          'stimulation_name', stimulation_name);
        ephysData(end+1) = newEntry;

        % Save the updated data to the MAT file
        save(matFilePath, 'ephysData');
        disp('Data added and saved.');
    else
        disp('This data combination already exists. No new data added.');
    end
end
