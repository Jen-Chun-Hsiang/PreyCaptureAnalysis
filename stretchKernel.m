function t_axis_stretched = stretchKernel(t_axis, kernel, stretch_factor)
    % Interpolates the kernel by a given stretch factor
    t_axis_stretched = linspace(t_axis(1)* stretch_factor, t_axis(end)* stretch_factor, length(kernel));
    %stretched_kernel = interp1(t_axis, kernel, t_axis_stretched, 'pchip', 'extrap');
end