function angleDegrees = calculateAngle(bodyLocation, objectLocation, bodyOrientationQuat)
    % Normalize the quaternion (important to ensure accuracy)
    bodyOrientationQuat = bodyOrientationQuat / norm(bodyOrientationQuat);

    % Convert the quaternion to a rotation matrix
    bodyRotationMatrix = quat2rotm(bodyOrientationQuat);

    % Calculate the object direction vector (from body to object)
    objectDirectionVector = objectLocation - bodyLocation;
    

    % Determine the forward direction of the body (assuming x-axis is forward)
    % This depends on how the coordinate system is defined in your application
    forwardVector = bodyRotationMatrix * [1; 1; 0];

    % Normalize the vectors
    objectDirectionVector = objectDirectionVector / norm(objectDirectionVector);
    forwardVector = forwardVector / norm(forwardVector);

    % Calculate the angle using the dot product
    cosTheta = dot(forwardVector, objectDirectionVector);
    angleRadians = acos(cosTheta);

    % Convert the angle to degrees
    angleDegrees = rad2deg(angleRadians);
end
