function registered = LineTransformationF(params, fixed, moving)
    % Extract parameters
    scale = params(1);
    translationX = params(2);

    % Create transformation
    tformMatrix = [scale, 0, 0; 0, 1, 0; translationX, 0, 1];
    tform = affine2d(tformMatrix);

    % Transform the moving image
    registered = imwarp(moving, tform, 'OutputView', imref2d(size(fixed)));

    % Calculate the similarity metric (mean squared error here)
    % val = immse(registered, fixed);
    % rmids = sum(registered, 1) == 0;
    % registered(:, rmids) = nan;
    % fixed(:, rmids) = [];
    % val = 1/(corr(registered(:), fixed(:)).^2+eps);
end