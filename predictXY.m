function [beta_x, beta_y, Y_pred] = predictXY(data, n, m)
    % data: Input data matrix where the first n columns are predictors,
    % and the last two columns are the outputs x and y to predict.
    % n: Number of predictor variables
    % m: Number of history steps for the predictors

    % Extracting the size of the data
    numRows = size(data, 1);

    % Initialize the feature matrix
    if m > 0
        numFeatures = n * m; % Total number of features with historical data
        Features = zeros(numRows - m, numFeatures);

        % Build feature matrix with lagged historical data for each predictor
        for i = 1:m
            for j = 1:n
                Features(:, n*(i-1)+j) = data(m+1-i:numRows-i, j);
            end
        end

        % Adjust the size of Target accordingly
        Target = data(m+1:end, end-1:end); % Targets are the last two columns
    else
        % When m = 0, use the current values of predictor variables directly
        Features = data(:, 1:n);
        Target = data(:, end-1:end);
    end

    % Adding a column of ones for the intercept in the regression model
    Features = [ones(size(Features, 1), 1) Features];

    % Fit linear regression models for x and y
    beta_x = regress(Target(:, 1), Features);
    beta_y = regress(Target(:, 2), Features);

    % Predict x and y using the fitted models
    Y_pred = Features * [beta_x beta_y];

end
