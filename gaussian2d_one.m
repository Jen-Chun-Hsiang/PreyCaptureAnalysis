function z = gaussian2d_one(x, y, params)
    x_mean = params(1);
    y_mean = params(2);
    sigma_x = params(3);
    sigma_y = params(4);
    theta = params(5);

    % Rotate coordinates
    x_rot = (x - x_mean) * cos(theta) + (y - y_mean) * sin(theta);
    y_rot = -(x - x_mean) * sin(theta) + (y - y_mean) * cos(theta);

    % Compute Gaussian
    z = exp(-(x_rot.^2/(2*sigma_x^2) + y_rot.^2/(2*sigma_y^2))) ;
end