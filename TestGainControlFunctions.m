
x = 0:0.1:10;
time_constants = [0.1 0.2 0.4 0.8 1.6];
num_time_constant = length(time_constants);



colors = parula(num_time_constant);
figure; hold on
for i = 1:num_time_constant
    y = decay(time_constants(i), x);
    plot(x, y, 'Color', colors(i, :));
end

%%
sigmoid = @(x, x0, tau) 1-1 ./ (1 + exp(-(x - x0) / tau));
x = 0:0.01:1;
figure; hold on
for i = 1:num_time_constant
    y = sigmoid(x, 0.5, time_constants(i));
    plot(x, y, 'Color', colors(i, :));
end
%%
power_gain = @(x, x0, tau) 1./(1+max([zeros(1, length(x)); (x-x0)], [], 1).^tau);
x = -1:0.1:10;
tau_list = 1:4;
num_tau = length(tau_list);
colors = parula(num_tau);
figure; hold on
for i = 1:num_tau
    y = power_gain(x, 0, tau_list(i));
    plot(x, y, 'Color', colors(i, :));
end
%%
% [History length, Sampling rate, Decay tau, Sigmoid tau, Sigmoid shift]
gain_params = [10, 10, 0.2, 0.2, 0.5];
function output = decay_func(tau, t) 
    output = exp(-t/tau);
end
function output = sigmoid_func(x, x0, tau) 
    output = 1-1 ./ (1 + exp(-(x - x0) / tau));
end
function G_t = compute_gain(params, history, t)
    G_t = sigmoid_func(decay_func(params(2), t)*history(:), params(4), params(3));
end
function [y, gt] = gain_control_system(params, x)
    history_len = params(1);
    y = nan(1, length(x)-history_len);
    gt = nan(1, length(x)-history_len);
    time_vec = (0:history_len-1)/params(2);
    tt = 1;
    for t = history_len+1:length(x)
        gt(t) = compute_gain(params(2:end), x(t-history_len:t-1), time_vec);
        y(t) = x(t-1)*gt(t);
        tt = tt + 1;
    end
end

%%
clc
gain_params = [20, 10, 0.3, 0.01, 0.01];
num_point = 100;
x = zeros(1, num_point);
x(51:70) = 0.01;
[y, gt] = gain_control_system(gain_params, x);
figure;subplot(1, 2, 1)
plot(x);
subplot(1, 2, 2);
hold on
plot(gt, 'k');
yyaxis right
plot(num_point-length(y)+1:num_point, y);


%% use the existing data to play arround the game control function
% run @CompareModelPrediction8Responses
c_dia_id = 6;
% sim = squeeze(resp_ON1(c_dia_id, 151:end));
% exp = squeeze(DataG(c_dia_id, 1, 151:end-1));
sim = squeeze(resp_OFF1(c_dia_id, 151:end));
exp = squeeze(DataG(c_dia_id, 2, 1:end-151));
figure; 
plot(exp);
yyaxis right
plot(sim);

%%
test_dia_ids = 3:7;
num_test_dia = length(test_dia_ids);
sim = [];
exp = [];
for j = 1:2
    for i = 1:num_test_dia
        switch j
            case 1
                csim = squeeze(resp_ON1(test_dia_ids(i), 151:end));
                cexp = squeeze(DataG(test_dia_ids(i), 1, 151:end-1));
            case 2
                csim = squeeze(resp_OFF1(test_dia_ids(i), 151:end));
                cexp = squeeze(DataG(test_dia_ids(i), 2, 1:end-151));
        end
        sim = [sim; csim(:)];
        exp = [exp; cexp(:)];
    end
end

ct = (0:length(sim)-1)/Fz;
figure; 
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, sim);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
box off
%%
% [History length, Sampling rate, Decay tau, Decay scalar, Sigmoid tau, Sigmoid shift]
gain_params = [50, 100, 0.3, 1, 0.01, 0.01];
CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], sim)*w(5)+w(6);...
    zeros(1, length(sim))], [], 1)-exp').^2, 'omitnan');
% Sigmoid
[OptW,fval] = fmincon(CostF, [0.1   1     0.01      0.01  1     0], [], [], [], [],...
                             [1e-8  0     1e-8      -2     1e-3  -100],...
                             [1    100   2         2     1e7   100]);
% Power
% [OptW,fval] = fmincon(CostF, [0.1   1     4   0.01   1     0], [], [], [], [],...
%                              [1e-8  0     1     -4   1e-3  -100],...
%                              [1     100   10     4   1e7   100]);

y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], sim)*OptW(5)+OptW(6);...
    zeros(1, length(sim))], [], 1);
ct = (0:length(sim)-1)/Fz;
figure; hold on
plot(ct, exp, 'Color', 0.5*ones(1, 3));
plot(ct, y, 'b');
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('Firing rate (spike/s)');

%%
dr_id = 3;
q = 1; %contra_id
i = 3; % bar_width_id
j = 1; % speed_id
sim = [];
exp = [];
for q = 1:3
    for i = 1:3
        csim = squeeze(resp(q, i, j, :));
        cexp = squeeze(mean(Data(dr_id, q, i, j, :, 1:length(csim)), 5, 'omitnan'));
        snan_ids = find(isnan(csim));
        enan_ids = find(isnan(cexp));
        if ~isempty(enan_ids)
            if snan_ids(1) > enan_ids(1)
                csim(snan_ids(1):end) = [];
                cexp(snan_ids(1):end) = [];
            else
                csim(enan_ids(1):end) = [];
                cexp(enan_ids(1):end) = [];
            end
        end
        snan_ids = find(isnan(csim));
        if ~isempty(snan_ids)
            csim(snan_ids) = csim(snan_ids(1)-1);
        end
        
        sim = [sim; csim];
        exp = [exp; cexp];
    end
end
ct = (0:length(sim)-1)/Fz;
figure; 
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, sim);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
%%
clear decay_func
function output = decay_func_show(tau, scalar, t)
    output = scalar*exp(-t/tau);
end
t = 0:0.1:(gain_params(1)/gain_params(2));

function output = sigmoid_func_show(x, x0, tau)
    output = 1-1 ./ (1 + exp(-(x - x0) / tau));
end

x = -10:0.001:10;

figure; 
subplot(1, 2, 1);
y = decay_func_show(OptW(1), OptW(2), t);
plot(t, y);
box off
xlabel('Time (s)');
title('Feedback filter');
subplot(1, 2, 2);
y = sigmoid_func_show(x, OptW(4), OptW(3));
plot(x, y);
box off
xlabel('Contrast');
title('Gain-control fucntion');
%%

t = 0:0.1:(gain_params(1)/gain_params(2));

function output = power_gain_show(x, x0, tau)
    output = 1./(1+max([zeros(1, length(x)); (x-x0)], [], 1).^tau);
end

x = -10:0.001:10;

figure; 
subplot(1, 2, 1);
y = decay_func_show(OptW(1), OptW(2), t);
plot(t, y);
box off
xlabel('Time (s)');
title('Feedback filter');
subplot(1, 2, 2);
y = power_gain_show(x, OptW(4), OptW(3));
plot(x, y);
box off
xlabel('Contrast');
title('Gain-control fucntion');