trajectory_name = 'trajectory_202404251111.mat';
trajectoryfolder = './Simulation/Trajectories';
slow_fac = 4;

spot_pos = load(fullfile(trajectoryfolder, trajectory_name), 'num_dots', 'positions');
num_dots = spot_pos.num_dots;
positions = spot_pos.positions;

x = 1:num_dots;
new_x = linspace(1, num_dots, num_dots*slow_fac);
new_num_dots = length(new_x);
new_positions = nan(new_num_dots, 2);

new_positions(:, 1) = interp1(x, positions(:, 1), new_x, 'pchip');
new_positions(:, 2) = interp1(x, positions(:, 2), new_x, 'pchip');

%%
Colors = parula(num_dots);
New_Colors = parula(new_num_dots);
close all

figure; hold on
for i = 1:(num_dots-1)
    plot(positions([i i+1], 1), positions([i i+1], 2), 'Color', Colors(i, :));
end

for i = 1:(new_num_dots-1)
    plot(new_positions([i i+1], 1), new_positions([i i+1], 2), '--', 'Color', New_Colors(i, :));
end
%%
keyboard;
%%
new_trajectory_name = sprintf('trajectory_202404251111_s%d.mat', slow_fac);
num_dots = new_num_dots;
positions = new_positions;
save(fullfile(trajectoryfolder, new_trajectory_name), 'num_dots', 'positions', 'slow_fac');
clear num_dots positions new_num_dots new_positions Colors New_Colors slow_fac spot_pos


