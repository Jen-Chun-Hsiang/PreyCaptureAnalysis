function [OptimizedParams, scaling] = GaussianTemporalFilter2(TargetFilters)
    y = TargetFilters;
    x = 1:length(y);
    scaling = max(abs(y));
    yb = y./scaling;

    [~, w1int] = max(yb);
    [~, w2int] = min(yb);

    % if w1int >= w2int
    %     cell_type = 'ON';
    % else
    %     cell_type = 'OFF';
    %     yb = -yb;
    % end


    % w = [sigma1, sigma2, mu1, mu2, amplitude, ratio, offset]
    CostF = @(w) mean( ...
        (gaussmf(x, [w(1) w(3)])*w(5) - gaussmf(x, [w(2) w(4)])*(w(6)*w(5)) + w(7) - yb).^2 ...
    );

    % Initial guess: [sigma1 sigma2 mu1 mu2 amplitude ratio offset]
    init = [2 2 w1int w2int -1 1 0];
    lb   = [0 0 min(x) min(x) -10 0  -1];
    ub   = [20 20 max(x) max(x)  10 10   1];

    [OptW, ~] = fmincon(CostF, init, [], [], [], [], lb, ub);

    % if strcmp(cell_type, 'OFF')
    %     OptW([5 7]) = -OptW([5 7]);
    % end
    OptimizedParams = OptW;
end