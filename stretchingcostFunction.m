function cost = stretchingcostFunction(stretch_factor, source_kernel, target_kernel, target_axis)
    source_axis = stretchKernel(target_axis, source_kernel, stretch_factor);
    left_t = max([source_axis(1) target_axis(1)]);
    right_t = min([source_axis(end) target_axis(end)]);
    new_axis = linspace(left_t, right_t, length(source_kernel));
    new_source = interp1(source_axis, source_kernel, new_axis, 'pchip', 'extrap');
    new_target = interp1(target_axis, target_kernel, new_axis, 'pchip', 'extrap');
    r = corr(new_source(:), new_target(:));
    cost = -r;
end