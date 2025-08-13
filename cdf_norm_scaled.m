function y = cdf_norm_scaled(x, mu, sigma, A, B)
  % Scaled and shifted cumulative normal distribution function
  z = (x - B) / abs(sigma); % standardize x
  y = A * normcdf(z, mu, 1) + B; % scale and shift the CDF
end