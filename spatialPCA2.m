function [eigvals_sorted, eigvecs_sorted, principal_components] = spatialPCA2(cData, locations)

points = locations;

% Step 1: Reshape the data to (T, n*m)
reshaped_data = cData;
num_image = size(cData, 1);

% Step 2: Construct a Spatial Weighting Matrix (W)

% Calculate pairwise Euclidean distances between sensor positions
distances = pdist2(points, points);

% Define a spatial weighting function (e.g., Gaussian function)
mean_distance = mean(distances(:));
W = exp(-distances.^2 / (2 * mean_distance^2));
% W = eye(length(W));

% Step 3: Modify the Covariance Matrix to include spatial weights
% Center the data
X_centered = reshaped_data - mean(reshaped_data, 1);

% Compute the weighted covariance matrix
weighted_cov = (X_centered' * X_centered)*W / num_image;

% is_symmetric = issymmetric(weighted_cov);

% Step 4: Perform Eigen Decomposition to get the spatial principal components
[eigvecs, eigvals] = eig(weighted_cov);

% Sort eigenvalues and eigenvectors in descending order
[eigvals_sorted, idx] = sort(diag(eigvals), 'descend');
eigvecs_sorted = eigvecs(:, idx);
eigvecs_sorted = real(eigvecs_sorted);
eigvals_sorted = eigvals_sorted*100/sum(eigvals_sorted);
eigvals_sorted = real(eigvals_sorted);

% Get the principal components
principal_components = X_centered * eigvecs_sorted;
% principal_components = principal_components';

end

