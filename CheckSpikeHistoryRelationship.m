clear all; close all; clc;
%% specify the recording neurons
load_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingWhite\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
recording_name = 'c081224';
load_recording_name = [recording_name '01'];
Fz = 100;
%% get the linear and nonlinear parts
load([load_data_folder load_recording_name '.mat'], 'masked_STAmat', 'PBs', 'FRs', 'stdSTA');

%%
timelapse = 100;
data1 = nan(1, timelapse);
for i = 1:timelapse
    B = FRs(1:end-i);
    A = FRs((i+1):end);
    A_given_B = A(B==1);
    data1(i) = sum(A_given_B)/sum(B==1);
end
data1 = data1-median(data1);
figure;
plot(data1)
%%
timelapse = 100;
data2 = nan(1, timelapse);
for i = 1:timelapse
    B = PBs(1:end-i);
    A = PBs((i+1):end);
    data2(i) = corr(A, B);
end
figure;
plot(data2)
%%
figure; 
plot(data1./(data2+eps));


%% get the spike data first from any stimulus processing file
fs = 10000;
Fz = 1000;
signals = zeros(size(sectionData));
signals(locs) = 1;
% t = (0:length(sectionData)-1)/fs;
downsample_size = round(fs/Fz);
add_more = mod(length(signals), downsample_size);
if add_more ~= 0
    add_more = downsample_size - add_more;
    sig = [signals(:)' zeros(1, add_more)];
end
assert(mod(length(sig), downsample_size) == 0);
sig = reshape(sig(:), downsample_size, []);
sig = sum(sig, 1);

%%
ssig = smoothdata(sig, 'movmean', 20);
len_sig = length(ssig);
time_gap = 5;
timelapse = 1:time_gap:500;
n_time = length(timelapse);
uni_v = unique(ssig);
n_uni_v = length(uni_v);
y = nan(n_uni_v, n_time);
for i = 1:n_time
    for k = 1:time_gap
        bids = 1:(len_sig-timelapse(i));
        aids = (timelapse(i)+1):len_sig;
        bids = bids+k-1;
        aids = aids+k-1;
        rmids = bids>len_sig | aids>len_sig;
        bids(rmids) = [];
        aids(rmids) = [];
        B = ssig(bids);
        A = ssig(aids);


        for j = 1:n_uni_v
            y(j, i, k) = mean(A(B==uni_v(j)), 'omitnan');
        end
    end

    % figure(1);
    % scatter(uni_v, y, 5, 'k', 'filled');
    % title(sprintf('- %d ms', timelapse(i)));
    % keyboard;
    % 
end
%%
plot_id = [1 2 12 13 14];
nplot = length(plot_id);
colors = parula(nplot);
figure; hold on
for i = 1:nplot
    plot(timelapse, mean(y(plot_id(i), :, :), 3, 'omitnan')*Fz, 'Color', colors(i, :));
    leg{i} = sprintf('%d', uni_v(plot_id(i))*Fz);
end
legend(leg);
ylabel('Firing rate (spike/s)');
%%
timelapse = 500;
data1 = nan(1, timelapse);
for i = 1:timelapse
    B = sig(1:end-i);
    A = sig((i+1):end);
    A_given_B = A(B==1);
    data1(i) = sum(A_given_B)/sum(B==1);
end
data1 = data1./median(data1);
figure;
plot(smoothdata(data1, 'movmedian', 10))
xlabel('ms');
ylabel('Probability');
