function [largest_segment, largest_area] = largest_segment_4conn_mask(binary_image)

% Validate input
if ~ismatrix(binary_image) || ~islogical(binary_image)
  error('Input must be a 2D logical matrix.');
end

% Find connected components and their properties
CC = bwconncomp(binary_image, 4); % 4-connected components
areas = cellfun(@numel, CC.PixelIdxList); % Area of each component

% Find the segment with the largest area
[largest_area, idx] = max(areas);
% Create a mask for the largest segment
largest_segment = zeros(size(binary_image));
largest_segment(CC.PixelIdxList{idx})= 1; % Pixels of largest segment

end
