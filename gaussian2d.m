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