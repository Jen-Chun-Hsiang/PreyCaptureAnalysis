function tf = gaussian_temporalfilter(x, OptW)
tf = (gaussmf(x, [OptW(1) OptW(3)])*OptW(5)-gaussmf(x, [OptW(2) OptW(4)])*OptW(6))+OptW(7);
