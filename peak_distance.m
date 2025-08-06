function [x, y, dist2d] = peak_distance(image)

% Check if input is valid
if ~ismatrix(image) || ~isnumeric(image)
  error('Input must be a numeric 2D matrix.');
end

% Find the peak intensity and its location
[peak_value, peak_idx] = max(image(:));
[peak_y, peak_x] = ind2sub(size(image), peak_idx);

% Calculate distance from peak for each pixel
[row_idx, col_idx] = ndgrid(1:size(image, 1), 1:size(image, 2));
distances = sqrt((row_idx - peak_y).^2 + (col_idx - peak_x).^2);

% Reflect distances across the peak
dist2d = reshape(distances, size(image));

% Combine distance, reflected distance, and pixel values
x = distances;
y = image;

end
