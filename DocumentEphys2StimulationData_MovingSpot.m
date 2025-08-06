clear; close all; clc
%%
ephys_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\ephys\sections\';
stim_data_folder = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Data\Stimulation\';
save_data_folder ='\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\Results\MovingSpot\';

for rr = 4 % reexamine 8 [6:7 12:40]
    recording_id = rr;
    switch recording_id
        case 1
            recordingname = 'a120723';
            load([ephys_data_folder recordingname '_0010_2.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingSpot_1.mat']);
            MinPeakHeight = 15;
        case 10
            recordingname = 'b120723';
            load([ephys_data_folder recordingname '_0002_2.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingSpot_1.mat']);
            MinPeakHeight = 20;
        case 11
            recordingname = 'c120723';
            load([ephys_data_folder recordingname '_0005_2.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingSpot_1.mat']);
            MinPeakHeight = 20;
        case 2
            recordingname = 'd010224';
            load([ephys_data_folder recordingname '_0001_2.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_2_MovingSpot_1.mat']);
        case 3
            recordingname = 'a030124';
            load([ephys_data_folder recordingname '_0000_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 90;
        case 4
            recordingname = 'd030124';
            load([ephys_data_folder recordingname '_0001_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 300;
        case 5
            recordingname = 'e030124';
            load([ephys_data_folder recordingname '_0000_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_000_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 200;
        case 6
            recordingname = 'a030924';
            load([ephys_data_folder recordingname '_0001_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_001_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 250;
        case 7
            recordingname = 'b030924';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 100;
        case 8
            recordingname = 'b030924';
            load([ephys_data_folder recordingname '_0003_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 100;
        case 9
            recordingname = 'c030924';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 200;
        case 39
            recordingname = 'd030924';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 400;
        case 12
            recordingname = 'a032924';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 250;
        case 13
            recordingname = 'b032924';
            load([ephys_data_folder recordingname '_0005_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_005_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 250;
        case 14
            recordingname = 'c032924';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 400;

        case 15 % not good detection
            recordingname = 'a040124';
            load([ephys_data_folder recordingname '_0003_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 500;
        case 16
            recordingname = 'b040124';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 500;
        case 34
            recordingname = 'c040124';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1500;
        case 35
            recordingname = 'd040124';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 2000;
        case 36
            recordingname = 'e040124';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 3000;
        case 37
            recordingname = 'f040124';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 2000;
        case 38
            recordingname = 'g040124';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1000;

        case 17  % bad detection
            recordingname = 'a040224';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 200;
        case 18
            recordingname = 'b040224';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 600;
        case 19
            recordingname = 'c040224';
            load([ephys_data_folder recordingname '_0005_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_005_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 300;
        case 20
            recordingname = 'd040224';
            load([ephys_data_folder recordingname '_0003_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 400;
        case 21
            recordingname = 'e040224';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 400;
        case 22
            recordingname = 'a040524';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 180;
        case 23
            recordingname = 'b040524';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1000;
        case 24
            recordingname = 'c040524';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 300;
        case 25
            recordingname = 'd040524';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 350;
        case 26
            recordingname = 'e040524';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 400;
        case 27
            recordingname = 'a040624';
            load([ephys_data_folder recordingname '_0003_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 700;
        case 28
            recordingname = 'b040624';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1000;
        case 29
            recordingname = 'c040624';
            load([ephys_data_folder recordingname '_0003_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1000;
        case 30
            recordingname = 'd040624';
            load([ephys_data_folder recordingname '_0003_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_003_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1000;
        case 31
            recordingname = 'e040624';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1000;
        case 32
            recordingname = 'f040624';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 500;
        case 33
            recordingname = 'g040624';
            load([ephys_data_folder recordingname '_0002_1.mat']);
            load([stim_data_folder 'Temporal_AlphaRGC_' recordingname '_002_Retina_1_MovingSpot.mat']);
            MinPeakHeight = 1000;
    end

    %%
    OUT = MS_OUT;
    IN = MS_IN;
    OLED.height = 600;
    OLED.width = 800;
    OLED.pixelSize = 3.5; %size of each pixel in micron on retina
    clear MS_OUT MS_IN
    assert(abs(length(sectionData)/(OUT.EndTime-OUT.StartTime) -10000)<1);


    %% Detrend
    close all
    figure;
    subplot(2, 1, 1); hold on
    plot(sectionData, 'k');
    % plot(smoothdata(sectionData, 'movmedian', 1*10000), 'r');
    lpass = lowpass(sectionData, 0.05, 10000, 'Steepness', 0.99);
    plot(lpass, 'r');

    subplot(2, 1, 2); hold on
    plot(sectionData, 'k');
    plot(sectionData-lpass, 'r');
    %%
    % keyboard;
    %%
    sectionData = sectionData-lpass;
    %%
    Fz = 100;
    WinT = [-0.5 0];
    % convert trace to spike or firing rate
    figure; hold on
    if ~exist('Is_extracellular', 'var')
        Is_extracellular = 1;
    end
    if Is_extracellular
        [pks, locs] = findpeaks(-sectionData, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 40);
        plot(sectionData(1:locs(end)), 'k');
        plot(locs(1:end), -pks(1:end), 'rx');
    else
        [pks, locs] = findpeaks(sectionData, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', 40);
        plot(sectionData(1:locs(1000)), 'k');
        plot(locs(1:1000), pks(1:1000), 'rx');
    end

    %%
    fs = 10000;
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
    t = (0:length(sig)-1)/Fz;
    ssigs = smoothdata(sig, 'gaussian', 0.05, 'SamplePoints', t');
    figure; hold on
    plot(t, sig, 'k');
    plot(t, ssigs, 'b');
    sig = ssigs*Fz;
    clear ssigs
    %%
    % keyboard;
    %% get new table with time ids
    % [Direction, Speed, Diameter]
    OUT.DataTable(:, 6) = round(Fz*(OUT.DataTable(:, 5)-OUT.StartTime))+1;
    Time_Difference_Check = ((OUT.EndTime-OUT.StartTime)*Fz-max(OUT.DataTable(:, 6)));
    % assert(((OUT.EndTime-OUT.StartTime)*Fz-max(OUT.DataTable(:, 6)))<=(IN.ISI*Fz+1));
    %%
    direc_ids = unique(OUT.DataTable(:, 4));
    speed_ids = unique(OUT.DataTable(:, 3));
    diame_ids = unique(OUT.DataTable(:, 2));
    Dtab = nan(length(direc_ids)*length(speed_ids)*length(diame_ids), 5);
    ii = 1;
    for i = 1:length(direc_ids)
        for j = 1:length(speed_ids)
            for k = 1:length(diame_ids)
                ctime = OUT.DataTable(OUT.DataTable(:, 4) == direc_ids(i) & OUT.DataTable(:, 3) == speed_ids(j) &...
                    OUT.DataTable(:, 2) == diame_ids(k) & ismember(OUT.DataTable(:, 1), [1 2]), 6);
                assert(length(ctime) == 2*IN.repeat);
                for v = 1:IN.repeat
                    Dtab(ii, :) = [i, j, k, ctime((v-1)*2+1), ctime(v*2)];
                    ii = ii + 1;
                end
            end
        end
    end
    Dtab(:, 6) = (Dtab(:, 5)-Dtab(:, 4))/Fz;
    %%
    MaxTwin = max(Dtab(:, 6))+IN.ISI; % in second (3.2 the max stimulus presentation length)

    %%
    addpath('\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\BasicFunctions');
    cData = nan(length(direc_ids), length(diame_ids), length(speed_ids));
    Data = nan(length(direc_ids), length(diame_ids), length(speed_ids), IN.repeat, round(MaxTwin*Fz));
    Moving_end_time_ids = nan(length(direc_ids), length(diame_ids), length(speed_ids), IN.repeat);
    clear Label_dias Label_spds
    close all
    for i = 1:length(direc_ids)
        figure;
        for k = 1:length(diame_ids)
            for j = 1:length(speed_ids)
                subplot(length(diame_ids), length(speed_ids), (k-1)*length(speed_ids)+j); hold on
                cids = find(Dtab(:, 1)==i & Dtab(:, 3)==k & Dtab(:, 2)==j);
                gsigs = nan(IN.repeat, round(MaxTwin*Fz));
                for v = 1:IN.repeat
                    ctim = Dtab(cids(v), 4:5);
                    csig = sig(ctim(1):ctim(2)+round(IN.ISI*Fz));
                    if length(csig) > round(MaxTwin*Fz)
                        gsigs(v, :) = csig(1:round(MaxTwin*Fz))';
                    else
                        gsigs(v, 1:length(csig)) = csig(:)';
                    end
                    ct = (0:length(csig)-1)/Fz;
                    plot(ct, csig, 'Color', 0.5*ones(1, 3));
                    cEndTid = ctim(2)-ctim(1)+1;
                    Moving_end_time_ids(i, k, j, v) = cEndTid;
                    plot(ct(cEndTid)*ones(1, 2), [0 60], 'r');
                    xlim([0, MaxTwin]);
                end
                Data(i, k, j, :, :) = gsigs;
                assert(sum(~isnan(gsigs(1, :))) > 1);
                cdata = gsigs(:, all(~isnan(gsigs), 1))';
                R2 = corr(cdata+eps*randn(size(cdata))).^2;
                cData(i, k, j) = meanR(R2);
                plot((0:size(gsigs, 2)-1)/Fz, mean(gsigs, 1), 'k');
                if k == 1
                    title(sprintf('Speeds:%d(um/s)', IN.speeds(j)))
                    Label_spds{j} = sprintf('%d (um/s)', IN.speeds(j));
                end
                if j == 1
                    if k ==1
                        ylabel(sprintf('Diameters:%d (um)', IN.diameters(k)));
                    else
                        ylabel(sprintf('%d (um)', IN.diameters(k)));
                    end
                    Label_dias{k} = sprintf('%d (um)', IN.diameters(k));
                elseif j == 2
                    ylabel('Firing rate (spike/s)')
                end
                ylim([0 150]);
            end
        end
        sgtitle(sprintf('%s Directions:%d', recordingname, IN.directions(i)));
    end
    %% Otherway to detect
    % (1) drop of the STD across repeat trials
    % (2) slope of accumulation of difference

    %% find the threshold for the response detection
    std_thr = 2.58; % prob (0.05 one tail)
    baseline_time = 0.3;
    signal_std = nan(length(direc_ids), 1);
    signal_mean = nan(length(direc_ids), 1);
    for i = 1:length(direc_ids)
        csignal = reshape(medfilt1(mean(Data(i, :, 1, :, 1:round(baseline_time*Fz)), 4), 5, [], 5), [], 1);
        signal_mean(i) = mean(csignal);
        signal_std(i) = std(csignal);
    end
    close all
    resp_onset = nan(length(direc_ids), length(diame_ids), length(speed_ids));
    t = (0:(size(Data, 5)-1))/Fz;
    for i = 1:length(direc_ids)
        figure;
        for k = 1:length(diame_ids)
            for j = 1:length(speed_ids)
                subplot(length(diame_ids), length(speed_ids), (k-1)*length(speed_ids)+j); hold on
                condition_signal =squeeze(Data(i, k, j, :, :));
                cal_type = 1;
                switch cal_type
                    case 1
                        avg_signal = medfilt1(mean(condition_signal, 1), 5);
                    case 2
                        avg_signal = std(medfilt1(condition_signal, 5, [], 2), [], 1);
                    case 3
                        avg_signal = medfilt1(mean(condition_signal, 1), 5);
                        avg_signal = cumsum(diff([avg_signal(1) avg_signal]));
                end
                % plot(t, condition_signal, 'Color', 0.5*ones(1, 3));
                plot(t, avg_signal, 'Color', 0*ones(1, 3));
                cids = find(abs(avg_signal-signal_mean(i)) > signal_std(i)*std_thr);

                stop_time = mean(Moving_end_time_ids(i, j, j, :)) ;
                if ~isempty(cids)
                    plot(t(cids(1)), avg_signal(cids(1)), 'xr');
                    resp_onset(i, k, j) = stop_time - cids(1);
                end
                plot(t(round(stop_time))*ones(1, 2), [0 100], 'b');
                xlim([0, MaxTwin]);
                if k == 1
                    title(sprintf('Speeds:%d(um/s)', IN.speeds(j)))
                    Label_spds{j} = sprintf('%d (um/s)', IN.speeds(j));
                end
                if j == 1
                    if k ==1
                        ylabel(sprintf('Diameters:%d (um)', IN.diameters(k)));
                    else
                        ylabel(sprintf('%d (um)', IN.diameters(k)));
                    end
                    Label_dias{k} = sprintf('%d (um)', IN.diameters(k));
                elseif j == 2
                    ylabel('Firing rate (spike/s)')
                end
                % ylim([0 150]);
            end
        end
        sgtitle(sprintf('%s Directions:%d', recordingname, IN.directions(i)));
    end

    %%
    % keyboard;
    %%
    dim1_moving_direction = IN.directions;
    dim2_diameters = IN.diameters;
    dim3_speeds = IN.speeds;
    dim4_repeats = 1:IN.repeat;
    dim5_time = (0:size(Data, 5)-1)/Fz;
    saveFileName = sprintf('%s_moving_spot_processed.mat', recordingname);
    save(sprintf('./Results/MovingSpot/%s', saveFileName), 'dim1_moving_direction', 'dim2_diameters', 'dim3_speeds',...
        'dim4_repeats', 'dim5_time', 'Data', 'Moving_end_time_ids', 'stop_time', 'std_thr', 'resp_onset', 'cData',...
        'recordingname', 'baseline_time','Time_Difference_Check');
    %%
    figure;
    imagesc(log(flipud(squeeze(mean(cData, 1))))); colorbar
    yticks(1:length(Label_dias));
    yticklabels(fliplr(Label_dias));
    xticks(1:length(Label_spds));
    xticklabels(Label_spds);
    title(sprintf('%s Response quality (log(R2)', recordingname));
end
