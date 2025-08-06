%%
dr_id = 3;
q = 1; %contra_id
i = 3; % bar_width_id
j = 1; % speed_id
ctrs = [1 2/3 1/3]*2;
sim = [];
exp = [];
ctr = [];
for k = 1:size(Data, 5)
    for q = 1:2
        iid = randperm(3);
        for i = 1:3
            csim = squeeze(resp(q, iid(i), j, :));
            cexp = squeeze(Data(dr_id, q, iid(i), j, k, 1:length(csim)));
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
            ctr = [ctr; ctrs(q)*ones(size(csim))];
        end
    end
end
ct = (0:length(sim)-1)/Fz;
contrast_name = 'makeup';
%% stimulus calculate of contrast
dr_id = 3;
q = 1; %contra_id
i = 3; % bar_width_id
j = 1; % speed_id
ctrs = [1 2/3 1/3]*2;
sim = [];
exp = [];
ctr = [];
for k = 1:size(Data, 5)
    for q = 1:1
        iid = randperm(3);
        for i = 1:3
            for j = 1:2
                csim = squeeze(resp(q, iid(i), j, :));
                ccntr = squeeze(cntr(q, iid(i), j, :));
                cexp = squeeze(Data(dr_id, q, iid(i), j, k, 1:length(csim)));
                snan_ids = find(isnan(csim));
                enan_ids = find(isnan(cexp));
                if ~isempty(enan_ids)
                    if snan_ids(1) > enan_ids(1)
                        csim(snan_ids(1):end) = [];
                        ccntr(snan_ids(1):end) = [];
                        cexp(snan_ids(1):end) = [];
                    else
                        csim(enan_ids(1):end) = [];
                        ccntr(enan_ids(1):end) = [];
                        cexp(enan_ids(1):end) = [];
                    end
                end
                snan_ids = find(isnan(csim));
                if ~isempty(snan_ids)
                    csim(snan_ids) = csim(snan_ids(1)-1);
                    ccntr(snan_ids) = ccntr(snan_ids(1)-1);
                end

                sim = [sim; csim];
                exp = [exp; cexp];
                ctr = [ctr; ccntr];
            end
        end
    end
end
ct = (0:length(sim)-1)/Fz;
contrast_name = 'std';
%%
figure; 
subplot(2, 1, 1)
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, sim);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
box off;

subplot(2, 1, 2); 
plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
hold on
plot(ct, ctr);
plot(ct, 1./ctr);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('Contrast status (arbi.)');
ylim([0 max(ctr(:))]);
box off;

%%
% [History length, Sampling rate, (1) Decay tau, (2) Decay scalar, (3) Sigmoid tau, (4) Sigmoid shift]
csim = sim;
csim = csim*1e6;
gain_params = [200, 100, 0.3, 1, 0.01, 0.01];
type_optimization = 'non';
is_sigmoid = 0;
switch type_optimization
    case 'non'
        if is_sigmoid
            CostF = @(w) mean((max([gain_control_system_opt_sigmoid([gain_params(1:2) w(1:4)], csim).*w(5)+...
                w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
        else
            CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*w(5)+...
                w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
        end
    case 'pos'
        CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*ctr*w(5)+...
            w(6).*ctr;zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
    case 'div'
        CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*(1./ctr)*w(5)+...
            w(6).*ctr;zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
    case 'div2'
        CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*(1./ctr)*w(5)+...
            w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
    case 'div3'
        CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*(1./ctr)*w(5)+...
            w(6).*(1./ctr);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
    case 'div4'
        CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*(1./(ctr+eps))*w(5)+...
            w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
    case 'div5'
        CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*(ctr+w(7))*w(5)+...
            w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
    case 'div6'
        CostF = @(w) mean((max([gain_control_system_opt([gain_params(1:2) w(1:4)], csim).*(1./(ctr+w(7)+eps))*w(5)+...
            w(6);zeros(1, length(csim))], [], 1)-exp').^2, 'omitnan');
end
rng('shuffle')
switch type_optimization
    
    case {'non', 'pos', 'div', 'div2', 'div3', 'div4'}
        if is_sigmoid
            [OptW,fval] = fmincon(CostF, ...
                [0.1   1    0.01   0.01  1     0    ], [], [], [], [],...
                [1e-8  0    1e-8   -2    1e-3  -200 ],...
                [1    1e3   2      2     1e7   200  ]);
        else
            [OptW,fval] = fmincon(CostF, ...
                [0.1   1    4   0.01    1     0    ], [], [], [], [],...
                [1e-8  0    1   -4      1e-3  -200 ],...
                [1    1e3   10   4      1e7   200  ]);
        end
    case {'div5', 'div6'}
        [OptW,fval] = fmincon(CostF,...
            [0.1   1    0.01   0.01  1     0    0.1], [], [], [], [],...
            [1e-8  0    1e-8   -2    1e-3  -200 0],...
            [1    1e3   2      2     1e7   200  2]);
end

switch type_optimization
    case 'non'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*OptW(5)+OptW(6);
            zeros(1, length(csim))], [], 1);
    case 'pos'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*ctr*OptW(5)+OptW(6).*ctr;
            zeros(1, length(csim))], [], 1);
    case 'div'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*(1./ctr)*OptW(5)+OptW(6).*ctr;
            zeros(1, length(csim))], [], 1);
    case 'div2'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*(1./ctr)*OptW(5)+OptW(6);
            zeros(1, length(csim))], [], 1);
    case 'div3'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*(1./ctr)*OptW(5)+OptW(6).*(1./ctr);
            zeros(1, length(csim))], [], 1);
    case 'div4'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*(1./(ctr+eps))*OptW(5)+OptW(6);
            zeros(1, length(csim))], [], 1);
    case 'div5'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*(ctr+OptW(7))*OptW(5)+OptW(6);
            zeros(1, length(csim))], [], 1);
    case 'div6'
        y = max([gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim).*(1./(ctr+OptW(7)+eps))*OptW(5)+OptW(6);
            zeros(1, length(csim))], [], 1);
end
ct = (0:length(sim)-1)/Fz;
figure; hold on
plot(ct, exp, 'Color', 0.5*ones(1, 3));
plot(ct, y, 'b');
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('Firing rate (spike/s)');
%%
foldername = './Results/Params';
save_file_name = sprintf('%s_%s_moving_bar_optimized_%d_opt-%s_contrast-%s_sigmoid%d_burry%d_%s.mat', recording_name,...
    response_name, implement_case_id, type_optimization, contrast_name, is_sigmoid, is_blurry, bar_type);
save(fullfile(foldername, save_file_name), 'OptW', 'Fz', 'y', 'sim', 'exp', 'type_optimization', 'gain_params');
%% Evaluate the parameters
figure; 
subplot(2, 1, 1)

plot(ct, exp);
ylabel('Firing rate (spike/s)');
yyaxis right
plot(ct, csim);
xlim([0 ct(end)]);
xlabel('Time (s)');
ylabel('linear RF signal (arbi.)');
box off;

subplot(2, 1, 2)
% yyaxis left
y = csim-min(csim(:));
y = y./range(y(:));
plot(ct, y);
ylabel('linear RF signal (arbi.)');
box off;
hold on

y = gain_control_system_opt([gain_params(1:2) OptW(1:4)], csim);
y = y-min(y(:));
y = y./range(y(:));
plot(ct, y);
ylabel('Gain control')

