blurry_length = 11;
figure; 
subplot(2, 3, 1);
imagesc(img'); colorbar
subplot(2, 3, 2);
bcr = max(img', -0.6*ones(size(img')));
maxv = max(bcr(:));
minv = min(bcr(:));
extend_ratio = (1-minv)/(maxv-minv);
bcr = bcr*extend_ratio+minv*(1-extend_ratio);
imagesc(bcr); colorbar
subplot(2, 3, 3);
bcrb = imgaussfilt(bcr, blurry_length);
imagesc(bcrb+img'); colorbar
subplot(2, 3, 4);
bcrba = bcrb+img';
maxv = max(bcrba(:));
minv = min(bcrba(:));
extend_ratio = (max(img(:))-min(img(:)))/(maxv-minv);
bcrba = bcrba*extend_ratio;
bcrba = bcrba-min(bcrba(:))+min(img(:));
imagesc(bcrba); colorbar
