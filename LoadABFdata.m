addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions\OpenSource');


sourceFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\';
destFolder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\';
convertAbfToMat(sourceFolder, destFolder);

function convertAbfToMat(sourceFolder, destFolder)
    % List all ABF files in the source folder
    abfFiles = dir(fullfile(sourceFolder, '*.abf'));

    for i = 1:length(abfFiles)
        abfFileName = abfFiles(i).name;
        matFileName = [abfFileName(1:end-4), '.mat']; % Change extension to .mat

        % Check if the MAT file already exists in the destination folder
        if ~exist(fullfile(destFolder, matFileName), 'file')
            % Load ABF file using abfload
            data = abfload(fullfile(sourceFolder, abfFileName));

            % Save the data to a MAT file in the destination folder
            save(fullfile(destFolder, matFileName), 'data');
            fprintf('Converted and saved: %s\n', matFileName);
        else
            fprintf('File already exists, skipped: %s\n', matFileName);
        end
    end
end
