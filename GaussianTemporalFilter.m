function [OptimizedParams, scaling] = GaussianTemporalFilter(TargetFilters)
    y = TargetFilters;
    x = 1:length(y);
    scaling = max(y);
    yb = y./scaling;


    [~, w1int] = max(yb);
    [~, w2int] = min(yb);
    CostF = @(w) mean((gaussian_temporalfilter(x, w) - yb).^2);
    [OptW,fval] = fmincon(CostF, [2    2  w1int   w2int     0  0  0], [], [], [], [],...
                                 [0    0  min(x)  min(x)    0  0 -1],...
                                 [20  20  max(x)  max(x)   10 10  1]);
    OptimizedParams = OptW;
    
    % yPred = (gaussmf(x, [OptW(1) OptW(3)])*OptW(5)-gaussmf(x, [OptW(2) OptW(4)])*OptW(6))+OptW(7);
    
%     figure;
%     T = (-20:1)./(9.47*2);
%     plot(T, yb, 'b'); hold on
%     plot(T, yPred, '-.b'); hold on
%     legend({' Noise',' Two Gaussian Fit'});
end