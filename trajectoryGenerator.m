function [positions, dot_poss] = trajectoryGenerator(num_pos, square_length, mass)
if nargin < 3
    mass = 10;
end

pos = tilingNeuronGenerator(num_pos, square_length);

velocities = zeros(num_pos, 2);
positions = zeros(num_pos, 2);
positions(1, :) = pos(1, :);
for i = 1:size(pos, 1)-1
    current_force = pos(i+1, :)-positions(i, :);
    dv = current_force/mass;
    velocities(i+1, :) = velocities(i, :) + dv;
    positions(i+1, :) = positions(i, :) + velocities(i+1, :);
end
% Centering
positions(:, 1) = positions(:, 1)*square_length/range(positions(:, 1));
positions(:, 2) = positions(:, 2)*square_length/range(positions(:, 2));
positions = positions - 0.5*(max(positions, [], 1)+min(positions, [], 1));
dot_poss = pos;

