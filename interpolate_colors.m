function color_list = interpolate_colors(color1, color2, num_colors)
%INTERPOLATE_COLORS Interpolates between two RGB colors.
%
% Inputs:
%   color1:     The start color (RGB vector, values between 0 and 1)
%   color2:     The end color (RGB vector, values between 0 and 1)
%   num_colors: The number of colors to generate (including start/end)
%
% Outputs:
%   color_list: A num_colors-by-3 matrix of interpolated RGB colors 

% Ensure colors are within the valid range (0 to 1)
if any(color1 < 0 | color1 > 1) || any(color2 < 0 | color2 > 1)
    error('Color values must be between 0 and 1');
end

% Create a vector of linearly spaced values between 0 and 1
interp_values = linspace(0, 1, num_colors);

% Interpolate each color channel independently
color_list = zeros(num_colors, 3);
for i = 1:3
    color_list(:, i) = interp1([0 1], [color1(i) color2(i)], interp_values);
end

end
