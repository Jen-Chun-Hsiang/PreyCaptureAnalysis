function connectivity = distanceConnection(target_pos, source_pos, rf_radius)
distance = sqrt((target_pos(:, 1)-source_pos(:, 1)').^2 + (target_pos(:, 2)-source_pos(:, 2)').^2);
connectivity = zeros(size(distance));
weight_fun = @(d, r) exp(-log(2)*d/(r/2));
for i = 1:size(distance, 1)
    cids = distance(i, :) < rf_radius;
    connectivity(i, cids) = weight_fun(distance(i, cids), rf_radius);
end
