function speeds = calculateSpeeds(xCoords, yCoords, frameRate)
    % Validate the input lengths
    if length(xCoords) ~= length(yCoords)
        error('xCoords and yCoords must be of the same length');
    end

    % Calculate differences between consecutive coordinates
    deltaX = diff(xCoords);
    deltaY = diff(yCoords);

    % Euclidean distances (displacements) between consecutive points
    displacements = sqrt(deltaX.^2 + deltaY.^2);

    % Time interval between frames
    timeInterval = 1 / frameRate; % seconds per frame

    % Speeds calculation
    speeds = displacements / timeInterval;
end