cell_id = 2;
% c = gaussian2d(X, Y, optimal_params(:, cell_id))';
c = double(MapMask == cell_id).*stdSTA;
c1 = c(:);
[c2, cids] = sort(c1, 'ascend');
d = cumsum(c2);
[~, rids] = sort(cids, 'ascend');
c1 = d(rids);
c1 = reshape(c1, size(c, 1), size(c, 2));

figure; 
subplot(1, 3, 1)
imagesc(c); colorbar
subplot(1, 3, 2)
imagesc(c1); colorbar
subplot(1, 3, 3)
imagesc(c1>(d(end)*0.5));colorbar

%%

