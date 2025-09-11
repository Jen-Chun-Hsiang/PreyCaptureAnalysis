function [sim_nl, LN_params] = fit_linear_transform_with_surround(sim, exp, sim_s, varargin)
    % FIT_LINEAR_TRANSFORM_WITH_SURROUND Fits a linear transformation with surround suppression and rectification
    %
    % Inputs:
    %   sim - simulated center data vector
    %   exp - experimental data vector
    %   sim_s - simulated surround data vector
    %   varargin - optional name-value pairs:
    %     'CSR_value' - target CSR (Center-Surround Ratio) value for gamma parameter
    %     'CSRStrength' - penalty weight for CSR constraint (default: 0, no penalty)
    %
    % Outputs:
    %   sim_nl - transformed simulation data (max(alpha*(sim - gamma*sim_s) + beta, 0))
    %   LN_params - structure containing fitted parameters
    %     .alpha - scaling factor for center-surround difference
    %     .beta  - offset
    %     .gamma - surround suppression strength (positive)
    
    % Parse optional inputs
    p = inputParser;
    addParameter(p, 'CSR_value', [], @isnumeric);
    addParameter(p, 'CSRStrength', 0, @isnumeric);
    parse(p, varargin{:});
    
    CSR_value = p.Results.CSR_value;
    CSRStrength = p.Results.CSRStrength;
    
    % Ensure inputs are column vectors
    sim = sim(:);
    exp = exp(:);
    sim_s = sim_s(:);
    
    % Remove any NaN values
    valid_idx = ~isnan(sim) & ~isnan(exp) & ~isnan(sim_s);
    sim_clean = sim(valid_idx);
    exp_clean = exp(valid_idx);
    sim_s_clean = sim_s(valid_idx);
    
    % Define cost function with center-surround and rectification nonlinearity
    % Include CSR penalty term if CSR_value and CSRStrength are provided
    if ~isempty(CSR_value) && CSRStrength > 0
        cost_func = @(params) mean((max(params(1) * (sim_clean - params(3) * sim_s_clean) + params(2), 0) - exp_clean).^2) + ...
                              CSRStrength * (params(3) + CSR_value)^2;
    else
        cost_func = @(params) mean((max(params(1) * (sim_clean - params(3) * sim_s_clean) + params(2), 0) - exp_clean).^2);
    end
    
    % Initial guess: alpha=1, beta=0, gamma=0.5 (or CSR_value if provided)
    if ~isempty(CSR_value)
        initial_gamma = CSR_value;
    else
        initial_gamma = 0.5;
    end
    initial_params = [1, 0, initial_gamma];
    
    % Set bounds to ensure gamma is positive
    lb = [-Inf, -Inf, 0]; % Lower bounds: gamma >= 0
    ub = [Inf, Inf, Inf];  % Upper bounds
    
    % Optimize parameters with bounds
    options = optimoptions('fmincon', 'Display', 'off');
    optimal_params = fmincon(cost_func, initial_params, [], [], [], [], lb, ub, [], options);
    
    % Extract fitted parameters
    LN_params.alpha = optimal_params(1);
    LN_params.beta = optimal_params(2);
    LN_params.gamma = optimal_params(3);
    
    % Store CSR penalty information
    LN_params.CSR_value = CSR_value;
    LN_params.CSRStrength = CSRStrength;
    if ~isempty(CSR_value)
        LN_params.CSR_deviation = abs(LN_params.gamma - CSR_value);
    else
        LN_params.CSR_deviation = NaN;
    end
    
    % Generate transformed simulation data with center-surround and rectification
    sim_nl = max(LN_params.alpha * (sim - LN_params.gamma * sim_s) + LN_params.beta, 0);
end
