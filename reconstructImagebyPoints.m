function rimg = reconstructImagebyPoints(img, x, y, v)
w = size(img, 1);
h = size(img, 2);

% Create grid for the image
[X, Y] = meshgrid(-w*0.5:(w*0.5-1), -h*0.5:(h*0.5-1));

rimg = griddata(x, y, v, X, Y, 'linear');
rimg_nearest = griddata(x, y, v, X, Y, 'nearest');
replace_ids = find(isnan(rimg));
[x_r, y_r] = ind2sub(size(rimg), replace_ids);
seed_ids = find(~isnan(rimg));
[x_s, y_s] = ind2sub(size(rimg), seed_ids);
rimg_nearest = griddata(x_s, y_s, rimg(seed_ids), x_r, y_r, 'nearest');
rimg(replace_ids) = rimg_nearest;
