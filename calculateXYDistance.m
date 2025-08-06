function xyDistance = calculateXYDistance(bodyLocation, objectLocation)
    % Extract the x and y coordinates
    bodyXY = bodyLocation(:, 1:2);
    objectXY = objectLocation(:, 1:2);

    % Calculate the distance in the x-y plane
    xyDistance = sqrt(sum((bodyXY - objectXY).^2, 2));
end