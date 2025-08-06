
% Step 1: Reshape the data to (T, n*m)
reshaped_data = cData;

% Step 2: Construct a Spatial Weighting Matrix (W)

% Calculate pairwise Euclidean distances between sensor positions
distances = pdist2(points, points);

% Define a spatial weighting function (e.g., Gaussian function)
mean_distance = mean(distances(:));
% W = exp(-distances.^2 / (2 * mean_distance^2));
W = eye(length(W));

% Step 3: Modify the Covariance Matrix to include spatial weights
% Center the data
X_centered = reshaped_data - mean(reshaped_data, 1);

% Compute the weighted covariance matrix
weighted_cov = (X_centered' * X_centered)*W / num_image;

% Step 4: Perform Eigen Decomposition to get the spatial principal components
[eigvecs, eigvals] = eig(weighted_cov);

% Sort eigenvalues and eigenvectors in descending order
[eigvals_sorted, idx] = sort(diag(eigvals), 'descend');
eigvecs_sorted = eigvecs(:, idx);

% Get the principal components
principal_components = eigvecs_sorted' * X_centered';

%%
% Step 5: Reshape and Analyze the Principal Components
% Reshape the first principal component back to (n, m)
first_pc = reshape(principal_components(1, :), [n, m]);
second_pc = reshape(principal_components(2, :), [n, m]);

% Visualization
figure;
imagesc(first_pc);
colorbar;
title('First Spatial Principal Component');

figure;
imagesc(second_pc);
colorbar;
title('Second Spatial Principal Component');
