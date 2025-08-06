nCh = size(data, 2);
figure; 
for i = 1:nCh
    subplot(nCh, 1, i);
    plot(data(:, i), 'k');
end