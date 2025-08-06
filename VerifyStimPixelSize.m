% Set the folder to search
folderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Shared\Andrew-Emily\Data\Stimulation'; % <-- Replace with your folder path

% Get list of .mat files matching the pattern
fileList = dir(fullfile(folderPath, '*Temporal_AlphaRGC_*_MovingNoise_*.mat'));


% Initialize output struct array
fileInfo = struct('name', {}, 'pixelSize', {});


for k = 1:length(fileList)
    fileName = fileList(k).name;
    filePath = fullfile(folderPath, fileName);

    % Extract string between "Temporal_AlphaRGC_" and "_MovingNoise_"
    expr = 'Temporal_AlphaRGC_(.*?)_MovingNoise_';
    tokens = regexp(fileName, expr, 'tokens');
    if ~isempty(tokens)
        extractedName = tokens{1}{1};
    else
        extractedName = '';
    end

    % Load OLED variable and get pixelSize, skip file if load fails
    try
        S = load(filePath, 'OLED');
        if isfield(S, 'OLED') && isfield(S.OLED, 'pixelSize')
            pxSize = S.OLED.pixelSize;
        else
            pxSize = NaN;
        end
    catch
        % Skip file if unable to load
        continue;
    end

    % Store in struct array
    fileInfo(end+1) = struct('name', extractedName, 'pixelSize', pxSize);
end


% Display results
disp('Extracted Names:');
disp({fileInfo.name});
disp('Pixel Sizes:');
disp([fileInfo.pixelSize]);

% Example data_set definition
data_sets = {'e100724', 'f100724', 'a101224', 'b101224',  'c101224',  'd101224', 'e101224',...
             'b101424', 'c101424', 'd101424', 'e101424',  'a101624',  'b101624', 'd101624', 'e101624',...
             'b101924', 'c101924', 'd101924', 'e101924',  'b103124',  'e103124', 'a110424', 'b110424',...
             'c110424', 'd110424', 'e110424', 'f110424',  'g110424',  'a110924', 'b110924', 'c110924',...
             'a111224'};

% Print names where pixelSize is not 4.375 and matches any entry in data_set
disp('Names with pixelSize not equal to 4.375 and matching data_set:');
for i = 1:length(fileInfo)
    if fileInfo(i).pixelSize ~= 4.375
        for j = 1:length(data_sets)
            if contains(fileInfo(i).name, data_sets{j})
                fprintf('Found matching name: %s with pixelSize: %.3f\n', fileInfo(i).name, fileInfo(i).pixelSize);
                break;
            end
        end
    end
end