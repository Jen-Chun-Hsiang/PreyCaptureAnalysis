function [normalized_freq, power_per_bin] = getSpatialFreqPower(img)

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
% figure; histogram(radial_freq(:));

% Radial binning of the power spectrum (averaging power over rings of constant frequency)
num_bins = 100;  % Number of bins for averaging
freq_bins = linspace(0, max_freq, num_bins);
power_per_bin = zeros(size(freq_bins));

% Average power in each frequency bin (1D radial averaging)

for i = 1:num_bins-1
    mask = radial_freq >= freq_bins(i) & radial_freq < freq_bins(i+1);
    % figure(1);
    % imagesc(mask); 
    % pause(0.5)
    power_per_bin(i) = mean(power_spectrum(mask), 'omitnan');  % Average power in this bin
end
% power_per_bin = power_per_bin./sum(power_per_bin);

% Normalize frequency bins to [0, 1] (optional)
normalized_freq = freq_bins / max_freq;
normalized_freq = normalized_freq(1:end-1);
power_per_bin = power_per_bin(1:end-1);