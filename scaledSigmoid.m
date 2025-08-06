function y = scaledSigmoid(x, L, k, x0, y0)
  y = L ./ (1 + exp(-k * (x - x0))) + y0;
end
