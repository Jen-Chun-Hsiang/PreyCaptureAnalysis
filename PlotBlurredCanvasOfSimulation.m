figure; hold on; plot(canvas(:, 100), 'Color', [140 23 230]/255);
%%
plot(canvas(:, 100), 'Color', [230 140 23]/255);
xlabel('Distance (pixel)');
ylabel('Contrast value');
title(sprintf('gaussian smooth: %0.3G', blurry_length));