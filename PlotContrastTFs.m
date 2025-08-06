%%
contrast_counts = nan(numContrasts, 1);
for i = 1:numContrasts
    contrast_counts(i) = sum(OUT.contrastPerBlock(NoiBIds) == contrastVals(i));
end
%%
STAmat_ct = STAmat_ct./reshape(contrast_counts, 1, 1, 1, []);
%%
colors = lines(numContrasts);
clear contrast_tag
figure; hold on
for i = 1:numContrasts
    contrast_id = i;
    
    smtstdSTA = medfilt2(stdSTA);
    [x, y, dist2d] = peak_distance(smtstdSTA);
    minD = 100;
    ythr = quantile(y(x>minD), 0.99);
    binary_image = smtstdSTA>ythr;
    [largest_segment_mask, ~] = largest_segment_4conn_mask(binary_image);
    
    smtSTAmat = medfilt3(squeeze(STAmat_ct(:, :, :, contrast_id)));
    masked_stdSTA = largest_segment_mask'.*medfilt2(std(smtSTAmat, [], 3))';
    smtSTAmat = permute(smtSTAmat, [2, 1, 3]);
    tRF = reshape(smtSTAmat, [], size(STAmat, 3))'*masked_stdSTA(:);
    t = WinT(1):1/Fz:WinT(end);
    tRF = (tRF-tRF(end))/range(tRF);
    %
    plot(t(2:end), tRF, 'Color', colors(i, :));
    xlabel('Time (s)');
    ylabel('STA contrast');
    box off
    title(save_recording_name);
    contrast_tag{i} =  sprintf('%2G', contrastVals(i));
end
legend(contrast_tag);
