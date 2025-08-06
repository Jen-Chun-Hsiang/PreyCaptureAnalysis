function [horizontalAngle, verticalAngle] = calculatePOV(objectPos, eyePos, eyeGazeQuaternion)
    % Normalize the quaternion
    eyeGazeQuaternion = eyeGazeQuaternion / norm(eyeGazeQuaternion);
    
    % Convert the quaternion to a rotation matrix
    eyeRotationMatrix = quat2rotm(eyeGazeQuaternion);

    % Calculate the vector from the eye to the object
    eyeToObjectVec = objectPos - eyePos;
    % Check for degenerate case (object at the same position as the eye)
    if norm(eyeToObjectVec) == 0
        horizontalAngle = 0;
        verticalAngle = 0;
        return;
    end
    
    % Rotate the vector into the eye's coordinate system
    try
        transformedVec = eyeRotationMatrix \ eyeToObjectVec';
    catch 
        keyboard;
    end
    
    % Calculate horizontal and vertical angles
    horizontalAngle = atan2(transformedVec(2), transformedVec(1));
    verticalAngle = atan2(transformedVec(3), sqrt(transformedVec(1)^2 + transformedVec(2)^2));

    % Convert angles to degrees
    horizontalAngle = rad2deg(horizontalAngle);
    verticalAngle = rad2deg(verticalAngle);
end