% Define parameters
g = -1;          % Coupling gain
lambda = 5;     % Coupling length constant
kernelSize = 100; % Size of the kernel (choose based on desired range)
[dx, dy] = meshgrid(-kernelSize:kernelSize, -kernelSize:kernelSize);
d = sqrt(dx.^2 + dy.^2);

% Compute the kernel
A = exp(-d / lambda);
A = A./sum(A, 'all');
K = -g * A;
% Exclude the center element from the sum
% K_center = -g;  % At d = 0, K = -g * exp(0) = -g
% sum_neighbors = sum(K(:)) - K_center;
% Set the center element
% K(kernelSize+1, kernelSize+1) = 1 - sum_neighbors;
figure; imagesc(K); colorbar
% Adjust the center element
K(kernelSize+1, kernelSize+1) = g;


%%
img_file_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\MEA\McGillDataset\MatImages';
file_name = 'FileTable.mat';
load(fullfile(img_file_folder, file_name), 'FileTable');
num_image = size(FileTable, 1);
i = 2;
file_name = sprintf('ResizeImg%d_%d.mat', FileTable(i, 1), FileTable(i, 2));
load(fullfile(img_file_folder, file_name), 'img');
img = squeeze(mean(double(img), 3))/255;
img = 2*(img-0.5);
% Apply the filter
rimg = imfilter(img, K, 'same', 'replicate');
figure; 
subplot(3, 1, 1);
imagesc(img);colorbar
title('original image');
subplot(3, 1, 2);
imagesc(rimg); colorbar
title('coupling addition');
subplot(3, 1, 3);
imagesc(rimg+img); colorbar
title('processed image');
sgtitle(sprintf('g:%0.3G, lambda:%0.3G', g, lambda))
