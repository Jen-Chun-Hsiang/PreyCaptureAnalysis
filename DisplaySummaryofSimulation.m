compound_datafolder = './Simulation/Data/Compound';
a = load(fullfile(compound_datafolder, 'onlyOFFRGCs.mat'));
a1 = load(fullfile(compound_datafolder, 'ON_1conn_OFFRGCs.mat'));
a4 = load(fullfile(compound_datafolder, 'ON_4conn_OFFRGCs.mat'));
a2 = load(fullfile(compound_datafolder, 'ON_2conn_OFFRGCs.mat'));
test_spot_pos = a.test_spot_pos;

figure; hold on
xlim([-500 500])
ylim([-500 500])
xlabel('Width')
ylabel('Height')

Colors = linspace(0, 0.6, size(test_spot_pos, 1))'.*ones(size(test_spot_pos, 1), 3);
for i = 1:size(test_spot_pos, 1)-1
    plot(test_spot_pos([i i+1], 1), test_spot_pos([i i+1], 2), 'Color', Colors(i, :));
end

Colors = interpolate_colors([12 0 255]/255, [158 159 237]/255, size(test_spot_pos, 1));
for i = 1:size(test_spot_pos, 1)-1
    plot(a1.Y_pred([i i+1], 1), a1.Y_pred([i i+1], 2), 'Color', Colors(i, :));
end

Colors = interpolate_colors([245 160 5]/255, [233 203 148]/255, size(test_spot_pos, 1));
for i = 1:size(test_spot_pos, 1)-1
    plot(a4.Y_pred([i i+1], 1), a4.Y_pred([i i+1], 2), 'Color', Colors(i, :));
end

% Colors = interpolate_colors([30 142 5]/255, [151 220 135]/255, size(test_spot_pos, 1));
% for i = 1:size(test_spot_pos, 1)-1
%     plot(a2.Y_pred([i i+1], 1), a2.Y_pred([i i+1], 2), 'Color', Colors(i, :));
% end

%%
Colors = interpolate_colors([181 0 82]/255, [204 141 170]/255, size(test_spot_pos, 1));
for i = 1:size(test_spot_pos, 1)-1
    plot(a.Y_pred([i i+1], 1), a.Y_pred([i i+1], 2), 'Color', Colors(i, :));
end


%%
figure; 
subplot(1, 3, 1)
hold on
plot(a1.summary_data_numlatent(:, 1), 'Color',[12 0 255]/255);
plot(a2.summary_data_numlatent(:, 1), 'Color',[30 142 5]/255);
plot(a4.summary_data_numlatent(:, 1), 'Color',[245 160 5]/255);
legend({'zero-connected', '1-connected', '3-connected'});
xlim([1 num_latent]);
ylim([0 1]);
xlabel('Principal components');
ylabel('Prediction correlation coefficient');

subplot(1, 3, 2)
hold on
plot(a1.summary_data_numlatent(:, 2), 'Color',[12 0 255]/255);
plot(a2.summary_data_numlatent(:, 2), 'Color',[30 142 5]/255);
plot(a4.summary_data_numlatent(:, 2), 'Color',[245 160 5]/255);
legend({'zero-connected', '1-connected', '3-connected'})
xlim([1 num_latent]);
ylim([0 1]);
xlabel('Principal components');
ylabel('Prediction correlation coefficient');

subplot(1, 3, 3)
hold on
plot(a1.summary_data_numlatent(:, 3), 'Color',[12 0 255]/255);
plot(a2.summary_data_numlatent(:, 3), 'Color',[30 142 5]/255);
plot(a4.summary_data_numlatent(:, 3), 'Color',[245 160 5]/255);
legend({'zero-connected', '1-connected', '3-connected'})
xlim([1 num_latent]);
xlabel('Principal components');
ylabel('Distance from target (um)');





