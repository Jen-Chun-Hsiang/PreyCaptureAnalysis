function episodeLabels = segmentEpisodes(timeStamps, frameTimeThreshold)
    % Calculate differences between consecutive timestamps
    timeDiffs = diff(timeStamps);

    % Create a binary array where 1 indicates a time gap within the threshold
    binaryArray = [1, timeDiffs <= frameTimeThreshold];

    % Use bwconncomp to find connected components (episodes)
    cc = bwconncomp(binaryArray);

    % Generate episode labels
    episodeLabels = zeros(1, length(timeStamps));
    for i = 1:cc.NumObjects
        episodeLabels(cc.PixelIdxList{i}) = i;
    end
end
