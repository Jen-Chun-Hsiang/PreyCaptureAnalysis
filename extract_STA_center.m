function cube_STA = extract_STA_center(masked_STAmat, center, mask_r_pix)
ex_ids = round((-mask_r_pix:mask_r_pix)+center(1));
ey_ids = round((-mask_r_pix:mask_r_pix)+center(2));
cube_STA = masked_STAmat(ex_ids, ey_ids, :);
figure;
imagesc(std(cube_STA, [], 3)); colorbar
end

% estimate the center
% resample the STA to fit to the center
% Cut the center cube, given the specified 