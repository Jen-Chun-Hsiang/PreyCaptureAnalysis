function [weights, y_pred] = calculateWeights(X, y)
    % CalculateWeights - Calculates regression weights using the normal equation
    % X is an m by n matrix of predictors
    % y is an m by 1 vector of observations

    % Adding a column of ones to X to account for the intercept
    X = [ones(size(X, 1), 1) X];
    
    % Calculating the weights using the normal equation
    % weights = (X'X)^(-1)X'y
    weights = pinv(X' * X) * X' * y;

    % Calculating predictions y' using the weights
    y_pred = X * weights;
end
