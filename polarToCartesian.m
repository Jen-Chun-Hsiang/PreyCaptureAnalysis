function [x, y] = polarToCartesian(r, theta)
    % Convert angle from degrees to radians
    theta_radians = deg2rad(theta);

    % Calculate x and y coordinates
    x = r * cos(theta_radians);
    y = r * sin(theta_radians);
end