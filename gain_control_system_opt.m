function y = gain_control_system_opt(params, x)
    history_len = params(1);
    y = nan(1, length(x)-history_len);
    gt = nan(1, length(x)-history_len);
    time_vec = (0:history_len-1)/params(2);
    tt = 1;
    for t = history_len+1:length(x)
        gt(t) = compute_gain(params(3:end), flipud(x(t-history_len:t-1)), time_vec);
        y(t) = x(t-1)*gt(t);
        tt = tt + 1;
    end

    function G_t = compute_gain(params, history, t)
        % G_t = sigmoid_func(decay_func(params(1), params(2), t)*history(:), params(4), params(3));
        G_t = power_gain(decay_func(params(1), params(2), t)*history(:), params(4), params(3));
    end
    function output = sigmoid_func(x, x0, tau)
        output = 1-1 ./ (1 + exp(-(x - x0) / tau));
    end
    function output = power_gain(x, x0, tau)
        output = 1./(1+max([zeros(1, length(x)); (x-x0)], [], 1).^tau);
    end
    function output = decay_func(tau, scalar, t)
        output = scalar*exp(-t/tau);
    end
end
