function smoothedSignal = smoothTimeSeries(signal, samplingRate)
    % Window size in milliseconds
    windowSizeMs = 50;

    % Convert window size to number of samples
    windowSizeSamples = windowSizeMs * (samplingRate / 1000);

    % Standard deviation of the Gaussian kernel
    % Typically, the window size should cover ~6 standard deviations
    sigma = windowSizeSamples / 6;

    % Create a Gaussian filter kernel
    % The kernel should span enough samples to cover the Gaussian
    kernelSize = ceil(windowSizeSamples);
    if mod(kernelSize, 2) == 0
        kernelSize = kernelSize + 1; % Ensure kernel size is odd
    end
    range = -(kernelSize-1)/2:(kernelSize-1)/2;
    gaussianKernel = exp(-0.5 * (range / sigma).^2);
    gaussianKernel = gaussianKernel / sum(gaussianKernel); % Normalize

    % Apply the Gaussian filter
    smoothedSignal = conv(signal, gaussianKernel, 'same');
end
