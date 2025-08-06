function canvas = gaussianblubs(maxW, maxH, SpatialWavelength, pixelSize, seed_id, is_size_fixed)
if nargin == 4
    seed_id = round(rand(1)*10000);
    is_size_fixed = 0;
elseif nargin==5
    is_size_fixed = 0;
end
exd = 3;
stream = RandStream('mrg32k3a','seed',seed_id+1);
sigma = SpatialWavelength/(pixelSize*4);
fix_size = 50;
if is_size_fixed
    canvas = imgaussfilt(randn(stream, maxH+round(fix_size*exd*2), maxW+round(fix_size*exd*2)), sigma);
    exdlen = round(exd*fix_size);
else
    canvas = imgaussfilt(randn(stream, maxH+round(sigma*exd*2), maxW+round(sigma*exd*2)), sigma);
    exdlen = round(exd*sigma);
end

canvas = canvas((exdlen+1):(size(canvas, 1)-exdlen), (exdlen+1):(size(canvas, 2)-exdlen));
canvas = canvas-min(canvas(:));
canvas = canvas/max(canvas(:));