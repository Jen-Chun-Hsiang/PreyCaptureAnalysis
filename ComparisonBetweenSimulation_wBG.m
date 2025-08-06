clear; close all; clc;
%%
compound_datafolder = './Simulation/Data/Compound';
compound_save_file_name ='ONconnect4_0611_BG50.mat';
B50_S1_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0611_BG50.mat';
B50_S1_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect4_0607_BG50.mat';
B50_S2_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0607_BG50.mat';
B50_S2_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect4_0613_BG100.mat';
B100_S1_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0613_BG100.mat';
B100_S1_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect4_0620_BG100.mat';
B100_S2_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0620_BG100.mat';
B100_S2_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect4_0615_BG150.mat';
B150_S1_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0615_BG150.mat';
B150_S1_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect4_0618_BG150.mat';
B150_S2_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0618_BG150.mat';
B150_S2_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect4_0619_BG200.mat';
B200_S1_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0619_BG200.mat';
B200_S1_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect4_0621_BG200.mat';
B200_S2_D4 = load(fullfile(compound_datafolder, compound_save_file_name));
compound_save_file_name ='ONconnect1_0621_BG200.mat';
B200_S2_D1 = load(fullfile(compound_datafolder, compound_save_file_name));
pos_ON = B200_S1_D4.pos_ON;
pos_OFF = B200_S1_D4.pos_OFF;
if ~isfield(B200_S1_D1, 'pixel2um')
    pixel2um = 2.5;
else
    pixel2um = B200_S1_D1.pixel2um;
end

%%
Fz = 100;
figure; hold on 
t = (0:size(B50_S1_D1.bg_pos, 1)-1)/Fz;
plot(t, B50_S1_D1.bg_pos(:, 2), 'k');
plot(t, B50_S1_D1.spot_pos(:, 1), 'b');
ylabel('Position');
xlabel('Time (s)');
legend({'Background', 'Spot'});
%%
close all
%% Examine the responses of one cell
Dset = B200_S1_D4;
% Left side
% y = [61.1695 7.1249 -3.22175];
% x = [-91.2388 -193.623 -272.105];
% right side
y = [69.7529 -11.352 13.6076];
x = [36.2111 158.461 250.155];
normalize = @(x) (x-min(x))/range(x);
cids = nan(length(x), 1);
for i = 1:length(x)
    dists = sqrt((Dset.pos_NFs(:, 1)-x(i)).^2 + (Dset.pos_NFs(:, 2)-y(i)).^2);
    [~, sids] = sort(dists, 'ascend');
    cids(i) = sids(1);
end
Colors = [1, 0, 0; 
          0, 1, 0;
          0, 0, 1];

figure; 
subplot(2, 1, 1); hold on
scatter(Dset.pos_NFs(:, 1)*pixel2um, Dset.pos_NFs(:, 2)*pixel2um, 25, 'k', 'filled');
subplot(2, 1, 2); hold on
plot(t, Dset.bg_pos(:, 2)*pixel2um, 'y');
plot(t, Dset.spot_pos(:, 1)*pixel2um, 'k');
for i = 1:length(x)
    subplot(2, 1, 1); hold on
    scatter(Dset.pos_NFs(cids(i), 1)*pixel2um, Dset.pos_NFs(cids(i), 2)*pixel2um, 35, Colors(i, :));
    subplot(2, 1, 2); hold on
    a = Dset.sim_NF(cids(i), :);
    a = (normalize(a) -0.5)*240;
    a = a-a(t==0.5);
    plot(t, a, 'Color', Colors(i, :));
end
%%
normalize = @(x) (x-min(x, [], 2))./range(x, 2);
Dset = B150_S1_D1;
figure; 
subplot(2, 1, 1); hold on
plot(t, Dset.bg_pos(:, 2), 'y');
plot(t, Dset.spot_pos(:, 1), 'k');subplot(2, 1, 2); hold on

a = normalize(Dset.sim_NF);
[~, sids] = sort(Dset.pos_NFs(:, 1), 'ascend');
a = a(sids, :);
subplot(2, 1, 2); 
imagesc(a);
xlim([1 size(a, 2)]);
ylim([1 size(a, 1)]);

%% Simply look at the average position based on the firing rate
normalize = @(x) (x-min(x, [], 2))./range(x, 2);
Dset = B100_S1_D1;
Only_OFF = 0;
cids = abs(Dset.pos_NFs(:, 2)) < 120;
tids = t>0.5;
figure; 
% subplot(3, 1, 1);
% imagesc(Dset.canvas(:, size(Dset.canvas, 2)*0.5-400+1:size(Dset.canvas, 2)*0.5+400));colormap(gray);
subplot(2, 1, 1); hold on
plot(t(tids), Dset.bg_pos(tids, 2)*pixel2um, 'y');
plot(t(tids), Dset.spot_pos(tids, 1)*pixel2um, 'k');
% subplot(2, 1, 2); hold on
if ~Only_OFF
    a = normalize(Dset.sim_NF(cids, tids));
else
    a = Dset.RGC2NF_OFF*Dset.sim_traces_OFF;
    a = normalize(a(cids, tids));
end
[~, sids] = sort(Dset.pos_NFs(cids, 1), 'ascend');
a = a(sids, :);
b = Dset.pos_NFs(cids, 1);
b = b(sids);
b = b'*(a./sum(a, 1));
plot(t(tids), b*pixel2um, 'b');
yticks(-600:600:600);
yticklabels({'-600', '0', '600'});
subplot(2, 1, 2); 
imagesc(a);
xlim([1 size(a, 2)]);
ylim([1 size(a, 1)]);
yticks([1 20.5 40])
yticklabels({'-1000', '0', '1000'});
set(gca,'XTick',[])

sgtitle(sprintf('Number of connection: %d', Dset.num_connect-1));

%% Show responses traces across simulation size and repeats
close all
normalize = @(x) (x-min(x, [], 2))./range(x, 2);
spatial_wavelength = 200;
Dset1 = B200_S1_D1;
Dset4 = B200_S1_D4;
Only_OFF = 0;
cids = abs(Dset1.pos_NFs(:, 2)) < 120;
tids = t>0.5;
[~, sids] = sort(Dset.pos_NFs(cids, 1), 'ascend');
b = Dset1.pos_NFs(cids, 1);
Colors = [0 0 255;
          255 0 0;
          0 128 128]/255;
figure; 
subplot(2, 1, 1); hold on
plot(t(tids), Dset1.bg_pos(tids, 2)*pixel2um, 'y');
plot(t(tids), Dset1.spot_pos(tids, 1)*pixel2um, 'k');
a1 = normalize(Dset1.sim_NF(cids, tids));
a1 = a1(sids, :);
cb = b(sids);
cb = cb'*(a1./sum(a1, 1));
plot(t(tids), cb*pixel2um, 'Color', Colors(1, :));

a4 = normalize(Dset4.sim_NF(cids, tids));
a4 = a4(sids, :);
cb = b(sids);
cb = cb'*(a4./sum(a4, 1));
plot(t(tids), cb*pixel2um, 'Color', Colors(2, :));

a = Dset1.RGC2NF_OFF*Dset1.sim_traces_OFF;
a = normalize(a(cids, tids));
a = a(sids, :);
cb = b(sids);
cb = cb'*(a./sum(a, 1));
plot(t(tids), cb*pixel2um, 'Color', Colors(3, :));
title('Positions');
ylabel('x positions');

yticks(-600:600:600);
yticklabels({'-600', '0', '600'});

subplot(2, 1, 2); 
imagesc(a4-a1);
xlim([1 size(a, 2)]);
ylim([1 size(a, 1)]);
ylabel('x positions');
yticks([1 20.5 40])
yticklabels({'-1000', '0', '1000'});
xlabel('Time(s)');
set(gca,'XTick',[])
title('Response difference');
sgtitle(sprintf('Spatial wavelength: %d (um)', spatial_wavelength*2));

%% 
Dset = B150_S1_D1;
is_ON = 0;
Contrast = {'OFF', 'ON'};
if is_ON
    cids = abs(pos_ON(:, 2)) < 120;
else
    cids = abs(pos_OFF(:, 2)) < 120;
end
tids = t>0.5;
figure; 
subplot(2, 1, 1); hold on
plot(t(tids), Dset.bg_pos(tids, 2), 'y');
plot(t(tids), Dset.spot_pos(tids, 1), 'k');
% subplot(2, 1, 2); hold on
if is_ON
    a = normalize(Dset.sim_traces_ON(cids, tids));
    [~, sids] = sort(pos_ON(cids, 1), 'ascend');
    b = pos_ON(cids, 1);
else
    a = normalize(Dset.sim_traces_OFF(cids, tids));
    [~, sids] = sort(pos_OFF(cids, 1), 'ascend');
    b = pos_OFF(cids, 1);
end

a = a(sids, :);
b = b(sids);
b = b'*(a./sum(a, 1));
plot(t(tids), b, 'b');
subplot(2, 1, 2); 
imagesc(a);
xlim([1 size(a, 2)]);
ylim([1 size(a, 1)]);
sgtitle(sprintf('Contrast %s, Number of connection: %d', Contrast{is_ON+1}, Dset.num_connect-1));
%% Use only spot moving period to predict the location
t_period = [1.7 2.2];
cids = t>t_period(1) & t<t_period(2);
X = B50_S1_D4.sim_NF(:, cids)';
[coeff,score,latent,tsquared,explained,mu] = pca(X);
y = B50_S1_D4.spot_pos(cids, 1);
beta_x = regress(y, score(:, 1:3));
% X = B50_S1_D1.sim_NF(:, cids)';
y_pred = score(:, 1:3) * beta_x;
X = B50_S2_D4.sim_NF(:, cids)';
X = coeff(:, 1:3)'*(X-mu)';
y_pred_2 = X' * beta_x;
X = B50_S2_D4.sim_NF';
X = coeff(:, 1:3)'*(X-mu)';
y_pred_3 = X' * beta_x;
figure; hold on
plot(t(cids), y, 'k');
plot(t(cids), y_pred, 'b');
plot(t(cids), y_pred_2, 'm');
plot(t, y_pred_3, 'g');

%% Use only spot moving period to predict the location
t_period = [1.7 2.2];
cids = t>t_period(1) & t<t_period(2);
X = B50_S1_D1.sim_NF(:, cids)';
[coeff,score,latent,tsquared,explained,mu] = pca(X);
y = B50_S1_D1.spot_pos(cids, 1);
beta_x = regress(y, score(:, 1:3));
% X = B50_S1_D1.sim_NF(:, cids)';
y_pred = score(:, 1:3) * beta_x;
X = B150_S1_D1.sim_NF(:, cids)';
X = coeff(:, 1:3)'*(X-mu)';
y_pred_2 = X' * beta_x;
X = B150_S1_D1.sim_NF';
X = coeff(:, 1:3)'*(X-mu)';
y_pred_3 = X' * beta_x;
figure; hold on
plot(t(cids), y, 'k');
plot(t(cids), y_pred, 'b');
plot(t(cids), y_pred_2, 'm');
plot(t, y_pred_3, 'g');
%%
assert(sum(D4.pos_NFs - D1.pos_NFs, 'all')== 0);
%%
figure; 
subplot(1, 3, 1)
imagesc(D4.sim_NF); colorbar
subplot(1, 3, 2)
imagesc(D1.sim_NF); colorbar
subplot(1, 3, 3)
imagesc(D4.sim_NF-D1.sim_NF);colorbar