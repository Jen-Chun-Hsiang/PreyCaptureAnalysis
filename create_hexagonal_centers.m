function points = create_hexagonal_centers(xlim, ylim, target_num_centers, varargin)
    % Default parameter values
    defaultMaxIterations = 100;
    defaultNoiseLevel = 0.3;
    defaultRandSeed = [];
    defaultNumPositions = [];
    defaultPositionIndices = [];
    
    % Input parser setup
    p = inputParser;
    addParameter(p, 'MaxIterations', defaultMaxIterations, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'NoiseLevel', defaultNoiseLevel, @(x) isnumeric(x) && isscalar(x));
    addParameter(p, 'RandSeed', defaultRandSeed, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'NumPositions', defaultNumPositions, @(x) isempty(x) || (isnumeric(x) && isscalar(x)));
    addParameter(p, 'PositionIndices', defaultPositionIndices, @(x) isempty(x) || isnumeric(x));
    
    % Parse the input arguments
    parse(p, varargin{:});
    
    max_iterations = p.Results.MaxIterations;
    noise_level = p.Results.NoiseLevel;
    set_rand_seed = p.Results.RandSeed;
    num_positions = p.Results.NumPositions;
    position_indices = p.Results.PositionIndices;
    
    x_min = xlim(1);
    x_max = xlim(2);
    y_min = ylim(1);
    y_max = ylim(2);
    x_range = x_max - x_min;
    y_range = y_max - y_min;

    % Set the random seed if specified
    if ~isempty(set_rand_seed)
        rng(set_rand_seed);
    end

    % Estimate initial side length based on target number of centers
    approximate_area = x_range * y_range;
    approximate_cell_area = approximate_area / target_num_centers;
    side_length = sqrt(approximate_cell_area / (3 * sqrt(3) / 2));

    % Calculate horizontal and vertical spacing
    dx = side_length * sqrt(3);
    dy = side_length * 1.5;

    % Estimate number of columns and rows
    cols = ceil(x_range / dx);
    rows = ceil(y_range / dy);

    % Function to generate grid points with given spacing and offsets
    function points = generate_points_with_noise(dx, dy, offset_x, offset_y, noise_level)
        points = [];
        for row = 0:(rows-1)
            for col = 0:(cols-1)
                x = col * dx + x_min - offset_x;
                y = row * dy + y_min - offset_y;
                if mod(row, 2) == 1
                    x = x + dx / 2;  % Offset every other row
                end

                % Add noise
                x = x + (rand - 0.5) * 2 * noise_level * dx;
                y = y + (rand - 0.5) * 2 * noise_level * dy;

                if x >= x_min && x < x_max && y >= y_min && y < y_max
                    points = [points; x, y];
                end
            end
        end
    end

    % Calculate initial offsets to center the grid
    offset_x = (cols * dx - x_range) / 2;
    offset_y = (rows * dy - y_range) / 2;

    % Generate initial grid points
    points = generate_points_with_noise(dx, dy, offset_x, offset_y, noise_level);

    % Adjust grid spacing for a fixed number of iterations
    for i = 1:max_iterations
        if size(points, 1) > target_num_centers
            dx = dx * 1.01;
            dy = dy * 1.01;
        else
            dx = dx * 0.99;
            dy = dy * 0.99;
        end
        cols = ceil(x_range / dx);
        rows = ceil(y_range / dy);
        offset_x = (cols * dx - x_range) / 2;
        offset_y = (rows * dy - y_range) / 2;
        points = generate_points_with_noise(dx, dy, offset_x, offset_y, noise_level);
        if abs(size(points, 1) - target_num_centers) <= target_num_centers * 0.05  % 5% tolerance
            break;
        end
    end

    % If specific indices are provided, return only those positions
    if ~isempty(position_indices)
        points = points(position_indices(position_indices <= size(points, 1)), :);
    % If a specific number of positions is requested, randomly select that many positions
    elseif ~isempty(num_positions) && num_positions < size(points, 1)
        indices = randperm(size(points, 1), num_positions);
        points = points(indices, :);
    end
end
