function positions = tilingNeuronGenerator(num_neuron, square_length, initial_fac)
if nargin  < 3
    initial_fac = 4;
end
pos = rand(num_neuron*initial_fac, 2)*square_length;
[X, Y] = meshgrid(1:length(pos), 1:length(pos));
distance = sqrt((pos(:, 1)'-pos(:, 1)).^2 + (pos(:, 2)'-pos(:, 2)).^2);

up_ids = triu(ones(length(distance), length(distance)),1)==1;
distance = distance(up_ids);
X = X(up_ids);
Y = Y(up_ids);
ii = 1;
while length(unique(X)) > num_neuron & ii < num_neuron*(initial_fac-1)
    clc
    fprintf('generator progress... %d/%d', ii, num_neuron*(initial_fac-1))
    [~, sids] = sort(distance, 'ascend');
    rmids = X == X(sids(1));
    distance(rmids) = [];
    Y(rmids) = [];
    X(rmids) = [];
    ii = ii +1;
end

uids = unique(X);

positions = pos(uids, :);
% 
% 
% close all
% figure; 
% subplot(1, 2, 1)
% scatter(pos(:, 1), pos(:, 2), 5, 'k', 'filled');
% subplot(1, 2, 2)
% scatter(pos(uids, 1), pos(uids, 2), 5, 'k', 'filled');