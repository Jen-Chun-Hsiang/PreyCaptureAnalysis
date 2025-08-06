x = bin_PBs;
y = bin_FRs*Fz;
gprMdl = fitrgp(x, y);

%%
[y_pred, y_std] = predict(gprMdl, x);
figure; hold on
scatter(x, y, 5, 'k', 'filled');
scatter(x, y_pred, 5, 'b', 'filled');

%%
[y_pred, y_std] = predict(gprMdl, PBs);
figure; 
plot(y_pred, 'k');
