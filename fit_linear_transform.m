function [sim_nl, LN_params] = fit_linear_transform(sim, exp)
    % FIT_LINEAR_TRANSFORM Fits a simple linear transformation with rectification
    %
    % Inputs:
    %   sim - simulated data vector
    %   exp - experimental data vector
    %
    % Outputs:
    %   sim_nl - transformed simulation data (max(alpha*sim + beta, 0))
    %   LN_params - structure containing fitted parameters
    %     .alpha - scaling factor
    %     .beta  - offset
    
    % Ensure inputs are column vectors
    sim = sim(:);
    exp = exp(:);
    
    % Remove any NaN values
    valid_idx = ~isnan(sim) & ~isnan(exp);
    sim_clean = sim(valid_idx);
    exp_clean = exp(valid_idx);
    
    % Define cost function with rectification nonlinearity
    cost_func = @(params) mean((max(params(1) * sim_clean + params(2), 0) - exp_clean).^2);
    
    % Initial guess
    initial_params = [1, 0]; % alpha=1, beta=0
    
    % Optimize parameters
    options = optimset('Display', 'off');
    optimal_params = fminsearch(cost_func, initial_params, options);
    
    % Extract fitted parameters
    LN_params.alpha = optimal_params(1);
    LN_params.beta = optimal_params(2);
    
    % Generate transformed simulation data with rectification
    sim_nl = max(LN_params.alpha * sim + LN_params.beta, 0);
end