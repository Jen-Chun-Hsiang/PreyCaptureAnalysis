clear clc
M = readmatrix('./Eyetracking/trajectories/m002_scan.csv');
M = M(all(~isnan(M(:, 22:25)), 2), :);
%%
figure;
imagesc(M(:, 22:24)'); colorbar
%% Separate episode
episodeLabels = segmentEpisodes(M(:, 25)', 0.01);
eps_ids = unique(episodeLabels);
c = tabulate(episodeLabels);
rmeps_ids = c(c(:, 2) <100, 1);
eps_ids(eps_ids==0 | ismember(eps_ids, rmeps_ids)) = [];
neps = length(eps_ids);
%%
Colors = parula(7);
spdbins = linspace(0, 2000, 25);
angbins = linspace(0, 180, 30);
distbins = linspace(0, 50, 25);
%%
IsDisplay = 0; 
%%
samplingRate = 200;
framerate = 200;
for i = 1:neps
    close all
    cids = find(episodeLabels == eps_ids(i));
    % smooth the data of each episode
    cM = M(cids, :);
    cM(any(isnan(cM), 2), :) = [];
    if isempty(cM)
        continue
    end
    for j = 1:24
        cM(:, j) = smoothTimeSeries(cM(:, j), samplingRate);
        %         figure; hold on
        %         plot(M(cids, j), 'k');
        %         plot(cM(:, j), 'b');
        %         keyboard;
    end
    % remove the begining and end due to smooth anormaly
    cM = cM(6:end-5, :);
    ncid = size(cM, 1);
    
    % cricket location in the visual field (Left eye)
    horizontalAngle = nan(ncid, 1);
    verticalAngle = nan(ncid, 1);
    for j = 1:ncid
        objectPos = cM(j, 22:24);
        eyePos = cM(j, 1:3);
        eyeGazeQuaternion = cM(j, 4:7);
        [horizontalAngle(j), verticalAngle(j)] = calculatePOV(objectPos, eyePos, eyeGazeQuaternion);
        if mod(j, 100) == 0
            clc
            fprintf('progress...cricket to eye...(%d/%d) \n', j, ncid);
        end
    end
    visuallocs_h_l = horizontalAngle;
    visuallocs_v_l = verticalAngle;
    
    % cricket location in the visual field (Right eye)
    horizontalAngle = nan(ncid, 1);
    verticalAngle = nan(ncid, 1);
    for j = 1:ncid
        objectPos = cM(j, 22:24);
        eyePos = cM(j, 8:10);
        eyeGazeQuaternion = cM(j, 11:14);
        [horizontalAngle(j), verticalAngle(j)] = calculatePOV(objectPos, eyePos, eyeGazeQuaternion);
        if mod(j, 100) == 0
            clc
            fprintf('progress...cricket to eye...(%d/%d) \n', j, ncid);
        end
    end
    visuallocs_h_r = horizontalAngle;
    visuallocs_v_r = verticalAngle;
    
    % relative mouse head position related to the cricket
    horizontalAngle = nan(ncid, 1);
    verticalAngle = nan(ncid, 1);
    relativeAngle = nan(ncid, 1);
    for j = 1:ncid
        objectPos = cM(j, 22:24);
        headPos = cM(j, 15:17);
        headQuaternion = cM(j, 18:21);
        [horizontalAngle(j), verticalAngle(j)] = calculatePOV(objectPos, headPos, headQuaternion);
        relativeAngle(j) = calculateAngle(headPos, objectPos, headQuaternion);
        if mod(j, 100) == 0
            clc
            fprintf('progress...cricket to body...(%d/%d) \n', j, ncid);
        end
    end
    relativeloc = cM(:, 22:23)-cM(:, 15:16);
    distance = calculateXYDistance(cM(:, 15:17), cM(:, 22:24));
    cricketvel = calculateSpeeds(cM(:, 22), cM(:, 23), framerate);
    headlocs_h = horizontalAngle;
    headlocs_v = verticalAngle;
    % separate the approaching period to the idle period
    
    speeds_l = calculateSpeeds(visuallocs_h_l, visuallocs_v_l, framerate);
    speeds_l = [nan; speeds_l];
    speeds_r = calculateSpeeds(visuallocs_h_r, visuallocs_v_r, framerate);
    speeds_r = [nan; speeds_r];
    cricketvel = [nan; cricketvel];
    %%
    cids = relativeAngle<90 & distance<20;
    cricketSpdThr = 4;
    figure;
    % Left eye position
    subplot(2, 4, 1); hold on
    plot(visuallocs_h_l(~cids), visuallocs_v_l(~cids), 'b');
    plot(visuallocs_h_l(cids), visuallocs_v_l(cids), 'r');
    plot(visuallocs_h_l(cids & cricketvel > cricketSpdThr), visuallocs_v_l(cids & cricketvel > cricketSpdThr), 'g');
    xlabel('Horizontal (deg)');
    ylabel('Vertical (deg)');
    title('Left eye field');
    subplot(2, 4, 5); hold on
    h1 = histogram(speeds_l(~cids), spdbins);
    h1.Normalization = 'Probability';
    h1.EdgeColor = 'w';
    h1.FaceColor = 'b';
    h2 = histogram(speeds_l(cids), spdbins);
    h2.Normalization = 'Probability';
    h2.EdgeColor = 'w';
    h2.FaceColor = 'r';
    h3 = histogram(speeds_l(cids & cricketvel > cricketSpdThr), spdbins);
    h3.Normalization = 'Probability';
    h3.EdgeColor = 'w';
    h3.FaceColor = 'g';
    xlabel('Cricket moving speed (deg/s)');
    ylabel('Probability');
        
    % Right eye position
    subplot(2, 4, 2); hold on
    plot(visuallocs_h_r(~cids), visuallocs_v_r(~cids), 'b');
    plot(visuallocs_h_r(cids), visuallocs_v_r(cids), 'r');
    plot(visuallocs_h_r(cids & cricketvel > cricketSpdThr), visuallocs_v_r(cids & cricketvel > cricketSpdThr), 'r');
    title('Right eye field');
    xlabel('Horizontal (deg)');
    ylabel('Vertical (deg)');
    subplot(2, 4, 6); hold on
    h1 = histogram(speeds_r(~cids), spdbins);
    h1.Normalization = 'Probability';
    h1.EdgeColor = 'w';
    h1.FaceColor = 'b';
    h2 = histogram(speeds_r(cids), spdbins);
    h2.Normalization = 'Probability';
    h2.EdgeColor = 'w';
    h2.FaceColor = 'r';
    h3 = histogram(speeds_r(cids & cricketvel > cricketSpdThr), spdbins);
    h3.Normalization = 'Probability';
    h3.EdgeColor = 'w';
    h3.FaceColor = 'g';
    xlabel('Cricket moving speed (deg/s)');
    ylabel('Probability');
    
    % Body position
    subplot(2, 4, 3); hold on
    plot(relativeloc(~cids, 1), relativeloc(~cids, 2), 'b');
    plot(relativeloc(cids, 1), relativeloc(cids, 2), 'r');
    plot(relativeloc(cids & cricketvel > cricketSpdThr, 1), relativeloc(cids & cricketvel > cricketSpdThr, 2), 'g');
    plot(relativeloc(end, 1), relativeloc(end, 2), 'xk');
    xlabel('x (cm)');
    ylabel('y (cm)');
    title('Cricket relative location to head');
    subplot(2, 4, 7); hold on
    h1 = histogram(distance(~cids), distbins);
    %h1.Normalization = 'Probability';
    h1.EdgeColor = 'w';
    h1.FaceColor = 'b';
    h2 = histogram(distance(cids), distbins);
    %h2.Normalization = 'Probability';
    h2.EdgeColor = 'w';
    h2.FaceColor = 'r';
    h3 = histogram(distance(cids & cricketvel > cricketSpdThr), distbins);
    %h2.Normalization = 'Probability';
    h3.EdgeColor = 'w';
    h3.FaceColor = 'g';
    xlabel('Cricket distance (cm)');
    ylabel('Probability');
    
    % Body position
    subplot(2, 4, 4); hold on
    plot(headlocs_h(~cids), headlocs_v(~cids), 'b');
    plot(headlocs_h(cids), headlocs_v(cids), 'r');
    plot(headlocs_h(end), headlocs_v(end), 'xk');
    xlabel('Horizontal (deg)');
    ylabel('Vertical (deg)');
    title('Cricket direction to head');
    subplot(2, 4, 8); hold on
    h1 = histogram(relativeAngle(~cids), angbins);
    %h1.Normalization = 'Probability';
    h1.EdgeColor = 'w';
    h1.FaceColor = 'b';
    h2 = histogram(relativeAngle(cids), angbins);
    %h2.Normalization = 'Probability';
    h2.EdgeColor = 'w';
    h2.FaceColor = 'r';
    xlabel('Cricket direction (deg)');
    ylabel('Probability');
    %%
    keyboard;
end


