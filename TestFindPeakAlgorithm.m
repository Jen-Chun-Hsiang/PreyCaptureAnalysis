smtstdSTA = imgaussfilt(stdSTA,5);
smtstdSTA = medfilt2(smtstdSTA, [5, 5]);

% smtstdSTA = medfilt2(smtstdSTA, [15, 15]);
% smtstdSTA = imgaussfilt(smtstdSTA, 32);
close all
[x, y, dist2d] = peak_distance(smtstdSTA);

figure; 
subplot(2, 2, 1);
imagesc(stdSTA');
subplot(2, 2, 2);
imagesc(smtstdSTA');
subplot(2, 2, 3);
imagesc(dist2d');
subplot(2, 2, 4);
scatter(x, y, 5, 'k', 'filled');

figure; 
subplot(1, 2, 1);
a = smtstdSTA';
a = (a-min(a(:)))/range(a(:));
imshow(a);
subplot(1, 2, 2);
imagesc(smtstdSTA'); hold on
plot([400 400], [0 600], '--k');
plot([0 800], [300 300], '--k');

%
a = a*256;
thr = quantile(a(:), 0.99);
Cent = FastPeakFind(a, thr);
Cent = reshape(Cent, 2, [])';

for i = 1:size(Cent, 1)
    hold on
    scatter(Cent(:, 1), Cent(:, 2), 10, 'red', 'filled');
end


%%
% Load the 2D image
smtstdSTA = medfilt2(stdSTA, [3, 3]);
image = smtstdSTA';

mask = masked_stdSTA>=0;
%mask = ~isnan(masked_stdSTA);
% Define the objective function to minimize
objective_function = @(params) gaussian_difference(params, image, mask);

% Initial guess for parameters: mean and covariance matrix
% initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1]; % [x_mean, y_mean, sigma_x, sigma_y, theta]
initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1, 200, 200, 0.1];

% Minimize the objective function using fminsearch
optimal_params = fminsearch(objective_function, initial_params);

% Generate the Gaussian model using the optimal parameters
[X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
gaussian_model = gaussian2d(X, Y, optimal_params);

% Plot the original image and the fitted Gaussian model
figure;
subplot(1, 2, 1);
% imshow(image, []);
imagesc(image); colorbar
title('Original Image');

subplot(1, 2, 2);
imagesc(gaussian_model); colorbar
title('Fitted Gaussian Model');

R = corr(reshape(image(masked_stdSTA>0), [], 1), reshape(gaussian_model(masked_stdSTA>0), [], 1))^2;

% Define the Gaussian function
% function z = gaussian2d(x, y, params)
%     x_mean = params(1);
%     y_mean = params(2);
%     sigma_x = params(3);
%     sigma_y = params(4);
%     theta = params(5);
%     bias = params(6);
%     scale = params(7);
% 
%     % Rotate coordinates
%     x_rot = (x - x_mean) * cos(theta) + (y - y_mean) * sin(theta);
%     y_rot = -(x - x_mean) * sin(theta) + (y - y_mean) * cos(theta);
% 
%     % Compute Gaussian
%     z = scale*exp(-(x_rot.^2/(2*sigma_x^2) + y_rot.^2/(2*sigma_y^2))) + bias;
% end
function z = gaussian2d(x, y, params)
    x_mean = params(1);
    y_mean = params(2);
    sigma_x = params(3);
    sigma_y = params(4);
    theta = params(5);
    bias = params(6);
    scale = params(7);

    s_sigma_x = params(8);
    s_sigma_y = params(9);
    s_scale = params(10);

    % Rotate coordinates
    x_rot = (x - x_mean) * cos(theta) + (y - y_mean) * sin(theta);
    y_rot = -(x - x_mean) * sin(theta) + (y - y_mean) * cos(theta);

    % Compute Gaussian
    z = scale*exp(-(x_rot.^2/(2*sigma_x^2) + y_rot.^2/(2*sigma_y^2))) +...
        s_scale*exp(-(x_rot.^2/(2*s_sigma_x^2) + y_rot.^2/(2*s_sigma_y^2))) + bias;
end



% Define the objective function to minimize (difference between image and Gaussian model)
function error = gaussian_difference(params, image, mask)
    % Generate the Gaussian model
    [X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
    gaussian_model = gaussian2d(X, Y, params);
    % Compute the difference between the image and the Gaussian model
    % error = sum(sum((image - gaussian_model).^2));
    %error = sum((image(mask) - gaussian_model(mask)).^2);
    error = 1-corr(image(:), gaussian_model(:));
end
%%
%% Receptor field size
smtstdSTA = medfilt2(stdSTA, [3, 3]);
minD = 100;
[x, y, dist2d] = peak_distance(smtstdSTA);
ythr = quantile(y(x>minD), 0.99);
binary_image = smtstdSTA>ythr;
[largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
masked_stdSTA = largest_segment_mask'.*smtstdSTA';
%%
keyboard
%%
% Load the 2D image
smtstdSTA = medfilt2(stdSTA, [3, 3]);
image = smtstdSTA;

% Define the objective function to minimize
mask = masked_stdSTA>0;
objective_function = @(params) sum_gaussian_difference(params, image, mask);

% Initial guess for parameters: mean and covariance matrix
initial_params = [size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1]; % [x_mean, y_mean, sigma_x, sigma_y, theta]

% Minimize the objective function using fminsearch
optimal_params = fminsearch(objective_function, initial_params);

% Generate the Gaussian model using the optimal parameters
[X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
gaussian_model = gaussian2d(X, Y, optimal_params);

% Plot the original image and the fitted Gaussian model
figure;
subplot(1, 2, 1);
% imshow(image, []);
imagesc(image); colorbar
title('Original Image');

subplot(1, 2, 2);
imagesc(gaussian_model); colorbar
title('Fitted Gaussian Model');

R1 = corr(reshape(image(masked_stdSTA>0), [], 1), reshape(gaussian_model(masked_stdSTA>0), [], 1))^2;


% Define the objective function to minimize (difference between image and Gaussian model)
function error = sum_gaussian_difference(params, image, mask)
    % Generate the Gaussian model
    [X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
    gaussian_model = gaussian2d(X, Y, params);
    
    % Compute the difference between the image and the Gaussian model
    error = sum((image(mask) - gaussian_model(mask)).^2);
    %error = 1-corr(image(:), gaussian_model(:));
end
%%
% Load the 2D image
smtstdSTA = medfilt2(stdSTA, [3, 3]);
image = smtstdSTA;
num_gaussians = 3; % Change this to the desired number of Gaussians

% Define the objective function to minimize
objective_function = @(params) mixture_gaussian_difference(params, image, num_gaussians);

% Initial guess for parameters: mean and covariance matrix for each Gaussian

initial_params = [];
for i = 1:num_gaussians
    initial_params = [initial_params, size(image, 1)/2, size(image, 2)/2, 50, 50, 0, 0.1, 1]; % [x_mean, y_mean, sigma_x, sigma_y, theta]
end

% Minimize the objective function using fminsearch
optimal_params = fminsearch(objective_function, initial_params);

param_length = numel(initial_params)/num_gaussians;
% Generate the mixture Gaussian model using the optimal parameters
[X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
mixture_gaussian_model = zeros(size(image));
for i = 1:num_gaussians
    gaussian_model = gaussian2d(X, Y, optimal_params((i-1)*param_length+1:i*param_length));
    mixture_gaussian_model = mixture_gaussian_model + gaussian_model;
end

% Plot the original image and the fitted mixture Gaussian model
figure;
subplot(1, 2, 1);
imshow(image, []);
title('Original Image');

subplot(1, 2, 2);
imshow(mixture_gaussian_model, []);
title('Fitted Mixture Gaussian Model');



% Define the objective function to minimize (difference between image and mixture Gaussian model)
function error = mixture_gaussian_difference(params, image, num_gaussians)
    [X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
    param_length = numel(params)/num_gaussians;
    mixture_gaussian_model = zeros(size(image));
    for i = 1:num_gaussians
        gaussian_model = gaussian2d(X, Y, params((i-1)*param_length+1:i*param_length));
        mixture_gaussian_model = mixture_gaussian_model + gaussian_model;
    end
    %error = sum(sum((image - mixture_gaussian_model).^2));
    error = 1-corr(image(:), gaussian_model(:));
end
Data