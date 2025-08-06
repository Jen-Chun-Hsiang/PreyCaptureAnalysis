% Define the size of the image (height x width)
img_height = 100;
img_width = 100;

% Number of dots (known points)
num_dots = 10;

% Generate random positions (x, y) and gray values for the dots
x_dots = randi([1, img_width], num_dots, 1);  % Random x positions
y_dots = randi([1, img_height], num_dots, 1); % Random y positions
gray_values = randi([0, 255], num_dots, 1);   % Random gray values (0-255)

% Create grid for the image
[X, Y] = meshgrid(1:img_width, 1:img_height);

% Use griddata to interpolate the gray values from the known dots
% 'linear' interpolation for linear extrapolation
img_reconstructed = griddata(x_dots, y_dots, double(gray_values), X, Y, 'linear');

% Fill NaN values resulting from the interpolation (due to extrapolation)
img_reconstructed(isnan(img_reconstructed)) = 0;

% Display the reconstructed image
imshow(uint8(img_reconstructed));
title('Reconstructed Image using Linear Extrapolation from Dots');

%% Load an example image (or you can use your own image)
img = imread('cameraman.tif');  % Load image (2D array, grayscale)
img = double(img);  % Convert image to double for computation

% Compute the 2D Fourier transform of the image
F = fft2(img);

% Shift the zero frequency component to the center of the spectrum
F_shifted = fftshift(F);

% Compute the power spectrum (magnitude squared of the Fourier coefficients)
power_spectrum = abs(F_shifted).^2;

% Get the size of the image
[rows, cols] = size(img);

% Define frequency axes (in terms of normalized frequency)
u = (-cols/2:(cols/2-1))/cols;  % Frequency axis for x dimension
v = (-rows/2:(rows/2-1))/rows;  % Frequency axis for y dimension
[U, V] = meshgrid(u, v);

% Compute the radial frequency (distance from the origin in frequency space)
radial_freq = sqrt(U.^2 + V.^2);

% Maximum frequency is limited by the image dimensions
max_freq = max(radial_freq(:));

% Radial binning of the power spectrum (averaging power over rings of constant frequency)
num_bins = 100;  % Number of bins for averaging
freq_bins = linspace(0, max_freq, num_bins);
power_per_bin = zeros(size(freq_bins));
close all
% Average power in each frequency bin (1D radial averaging)
figure; 
imagesc(radial_freq);colorbar

for i = 1:num_bins-1
    mask = radial_freq >= freq_bins(i) & radial_freq < freq_bins(i+1);
    % figure(1);
    % imagesc(mask); 
    % pause(0.5)
    power_per_bin(i) = mean(power_spectrum(mask), 'omitnan');  % Average power in this bin
end

% Normalize frequency bins to [0, 1] (optional)
normalized_freq = freq_bins / max_freq;

% Plot the 1D power as a function of spatial frequency
figure;
plot(normalized_freq(1:end-1), log10(power_per_bin(1:end-1)), 'LineWidth', 2);
xlabel('Normalized Spatial Frequency');
ylabel('Power');
title('Power Spectrum as a Function of 1D Spatial Frequency');
grid on;

%%




