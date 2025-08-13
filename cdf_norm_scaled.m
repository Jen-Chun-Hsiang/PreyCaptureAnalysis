function y = cdf_norm_scaled(x, A, mu, sigma, offset)
    y = A * normcdf(x, mu, sigma) + offset;
end