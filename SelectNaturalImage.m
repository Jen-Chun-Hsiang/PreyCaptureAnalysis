clear; close all; clc;
img_file_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\MEA\McGillDataset\MatImages';
file_name = 'FileTable.mat';
load(fullfile(img_file_folder, file_name), 'FileTable');
num_image = size(FileTable, 1);
usability = nan(num_image, 1);
save_folder = './Results/Params';
save_file_name = 'NaturalImageUsability.mat';
if exist(fullfile(save_folder, save_file_name), 'file')
    load(fullfile(save_folder, save_file_name))
    work_on_ids = find(isnan(usability));
    num_work_on = length(work_on_ids);
else
    num_work_on = num_image;
end
if num_work_on > 0
    for i = 1:num_work_on
        file_name = sprintf('ResizeImg%d_%d.mat', FileTable(work_on_ids(i), 1), FileTable(work_on_ids(i), 2));
        load(fullfile(img_file_folder, file_name), 'img');
        img = squeeze(mean(double(img), 3))/255;

        figure(1);
        imshow(img);
        clc
        user_input = input('Enter usability evaluation between 0 and 3:');  % Prompt user to enter 'Y' or 'N'
        % Check if the input is a valid number between 0 and 5
        if isnumeric(user_input) && user_input >= 0 && user_input <= 5
            usability(work_on_ids(i)) = user_input;  % Store the valid input
        else
            fprintf('Invalid input, please enter a number between 0 and 5.\n');
            i = i - 1;  % Decrease index to repeat the loop for valid input
        end
        save(fullfile(save_folder, save_file_name), 'usability');
    end
end