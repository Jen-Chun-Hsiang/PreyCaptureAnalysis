function smoothed = gaussian_smooth_1d(rgc_time, kernel_size, sampling_rate, sigma)
%GAUSSIAN_SMOOTH_1D  Smooth a 1-D signal with a Gaussian kernel
%
%  smoothed = gaussian_smooth_1d(rgc_time, kernel_size, sampling_rate, sigma)
%
%  Inputs:
%    rgc_time      – [T×1] or [1×T] vector of firing-rate or binned spike counts
%    kernel_size   – length of the Gaussian kernel in samples (must be odd)
%    sampling_rate – samples per second (e.g. 100 for 10 ms bins)
%    sigma         – Gaussian σ in seconds (so σ_samples = sigma*sampling_rate)
%
%  Output:
%    smoothed      – same size as rgc_time, smoothed via conv(...,'same')

    % ensure column vector
    rgc_time = rgc_time(:);

    % convert σ (sec) to samples
    sigma_samp = sigma * sampling_rate;

    % build a symmetric kernel of length kernel_size
    half_k = (kernel_size-1)/2;
    x = -half_k : half_k;               % sample offsets
    gauss = exp( - (x.^2) / (2*sigma_samp^2) );
    gauss = gauss / sum(gauss);         % normalize area = 1

    % convolve and keep same length
    smoothed = conv(rgc_time, gauss, 'same');
end
