function gaussian_model = gaussian_multi(params, image, num_gauss)
    % Generate the Gaussian model
    [X, Y] = meshgrid(1:size(image,2), 1:size(image,1));
    params = reshape(params(:), [], num_gauss)';
    
    for i = 1:num_gauss
        if i == 1
            gaussian_model = gaussian2d(X, Y, params(i, :));
        else
            gaussian_model = gaussian_model + gaussian2d(X, Y, params(i, :));
        end
    end
    % Compute the difference between the image and the Gaussian model
    % error = sum(sum((image - gaussian_model).^2));
    %error = sum((image(mask) - gaussian_model(mask)).^2);
    % error = 1-corr(image(:), gaussian_model(:));
end