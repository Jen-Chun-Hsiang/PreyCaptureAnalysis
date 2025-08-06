function [x_indices, y_indices] = sample_indices_2d(array_size, num_of_points, varargin)
    % Default parameter values
    defaultXBound = [1, array_size(1)];
    defaultYBound = [1, array_size(2)];
    defaultRandSeed = [];

    % Input parser setup
    p = inputParser;
    addParameter(p, 'XBound', defaultXBound, @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'YBound', defaultYBound, @(x) isnumeric(x) && numel(x) == 2);
    addParameter(p, 'RandSeed', defaultRandSeed, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));

    % Parse the input arguments
    parse(p, varargin{:});

    x_bound = p.Results.XBound;
    y_bound = p.Results.YBound;
    rand_seed = p.Results.RandSeed;

    % Set the random seed if specified
    if ~isempty(rand_seed)
        rng(rand_seed);
    end

    % Ensure that bounds are within the array size
    x_bound(1) = max(1, x_bound(1));
    x_bound(2) = min(array_size(1), x_bound(2));
    y_bound(1) = max(1, y_bound(1));
    y_bound(2) = min(array_size(2), y_bound(2));

    % Calculate the possible x and y indices
    possible_x = x_bound(1):x_bound(2);
    possible_y = y_bound(1):y_bound(2);

    % Ensure num_of_points does not exceed the total available points
    total_possible_points = length(possible_x) * length(possible_y);
    num_of_points = min(num_of_points, total_possible_points);

    % Randomly sample the x and y indices
    [x_grid, y_grid] = meshgrid(possible_x, possible_y);
    indices = randperm(total_possible_points, num_of_points);
    x_indices = x_grid(indices);
    y_indices = y_grid(indices);
end
