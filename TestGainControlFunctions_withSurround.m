%%
dr_id = 3;
q = 1; %contra_id
i = 3; % bar_width_id
j = 1; % speed_id
sim = [];
sim_s = [];
exp = [];
for q = 1:3
    for i = 1:3
        csim = squeeze(resp(q, i, j, :));
        csim_s = squeeze(resp_s(q, i, j, :));
        cexp = squeeze(mean(Data(dr_id, q, i, j, :, 1:length(csim)), 5, 'omitnan'));
        snan_ids = find(isnan(csim));
        enan_ids = find(isnan(cexp));
        if ~isempty(enan_ids)
            if snan_ids(1) > enan_ids(1)
                csim(snan_ids(1):end) = [];
                csim_s(snan_ids(1):end) = [];
                cexp(snan_ids(1):end) = [];
            else
                csim(enan_ids(1):end) = [];
                csim_s(enan_ids(1):end) = [];
                cexp(enan_ids(1):end) = [];
            end
        end
        snan_ids = find(isnan(csim));
        if ~isempty(snan_ids)
            csim(snan_ids) = csim(snan_ids(1)-1);
            csim_s(snan_ids) = csim_s(snan_ids(1)-1);
        end
        
        sim = [sim; csim];
        sim_s = [sim_s; csim_s];
        exp = [exp; cexp];
    end
end
ct = (0:length(sim)-1)/Fz;
figure; 
subplot(1, 2, 1)
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, sim);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
box off
title('excitation');

subplot(1, 2, 2)
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, sim_s);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
box off
title('inhibition');

%%
% [History length, Sampling rate, Decay tau, Decay scalar, Sigmoid tau, Sigmoid shift]
gain_params = [50, 100, 0.3, 1, 0.01, 0.01];
CostF = @(w) mean((max([(gain_control_system_opt([gain_params(1:2) w(1:4)], sim)-...
w(7)*(gain_control_system_opt([gain_params(1:2) w(1:4)], sim_s)-w(8)))*w(5)+w(6);
    zeros(1, length(sim))], [], 1)-exp').^2, 'omitnan');
% Sigmoid
[OptW,fval] = fmincon(CostF, [0.1   1    0.01   0.01  1     0     0.1  0], [], [], [], [],...
                             [1e-8  0    1e-8   -2    1e-3  -100  0    -1],...
                             [1    100   2      2     1e7   100   100  1]);
% Power
% [OptW,fval] = fmincon(CostF, [0.1   1     4   0.01   1     0], [], [], [], [],...
%                              [1e-8  0     1     -4   1e-3  -100],...
%                              [1     100   10     4   1e7   100]);


y = max([(gain_control_system_opt([gain_params(1:2) OptW(1:4)], sim)-...
    OptW(7)*(gain_control_system_opt([gain_params(1:2) OptW(1:4)], sim_s)-OptW(8)))*OptW(5)+OptW(6);
    zeros(1, length(sim))], [], 1);
ct = (0:length(sim)-1)/Fz;
figure; hold on
plot(ct, exp, 'Color', 0.5*ones(1, 3));
plot(ct, y, 'b');
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('Firing rate (spike/s)');