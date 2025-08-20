function tf = diffgauss_tf(x, params)
    tf = gaussmf(x, [params(1) params(3)])*params(5) - gaussmf(x, [params(2) params(4)])*(params(6)*params(5)) + params(7);
end