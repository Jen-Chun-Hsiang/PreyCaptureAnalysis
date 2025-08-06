image_folder = './Spectrum';
case_id = 4;
switch case_id
    case 1
        image_name = 'UV filter 0517/20240517_155632.jpg';
        A = imread(fullfile(image_folder, image_name));
        image_name = 'Green filter 0517/20240517_155610.jpg';
        B = imread(fullfile(image_folder, image_name));
    case 2
        image_name = 'Green filter 0529/20240529_170809.jpg';
        A = imread(fullfile(image_folder, image_name));
        image_name = 'UV filter 0529/20240529_170807.jpg';
        B = imread(fullfile(image_folder, image_name));
    case 3
        image_name = 'Green filter 0529/20240529_170833.jpg';
        A = imread(fullfile(image_folder, image_name));
        image_name = 'UV filter 0529/20240529_170831.jpg';
        B = imread(fullfile(image_folder, image_name));
    case 4
        image_name = 'Green filter 0529/20240529_170842.jpg';
        A = imread(fullfile(image_folder, image_name));
        image_name = 'UV filter 0529/20240529_170840.jpg';
        B = imread(fullfile(image_folder, image_name));
end
%%
norm = @(x) x./max(x(:), [], "all");
norm1 = @(x) x./max(x(1:3000, :, :), [], "all");
norm2 = @(x) x./max(x(1:2200, :, :), [], "all");
imagesc_display = 0; % 0: imshow, 1: imagesc with colorbar
is_normalized = 0;
figure;
for i = 1:2
    for j = 1:3
        subplot(2, 3, (i-1)*3+j);
        switch i
            case 1
                if ~imagesc_display
                    if is_normalized
                        imshow(norm1(double(A(:, :, j))));
                    else
                        imshow(double(A(:, :, j))/255);
                    end
                else
                    imagesc(norm1(double(A(:, :, j))), [0 1]);colorbar
                end
            case 2
                if ~imagesc_display
                    if is_normalized
                        imshow(norm1(double(B(:, :, j))));
                    else
                        imshow(double(B(:, :, j))/255);
                    end
                else
                    imagesc(norm1(double(B(:, :, j))), [0 1]);colorbar
                end
        end
    end
end

%%
figure;
for i = 1:2
    subplot(1, 2, i);
    switch i
        case 1
            if ~imagesc_display
                imshow(norm1(double(A(:, :, 2))));
            else
                imagesc(norm1(double(A(:, :, 2))), [0 1]);colorbar
            end
        case 2
            if ~imagesc_display
                imshow(norm1(double(B(:, :, 3))));
            else
                imagesc(norm1(double(B(:, :, 3))), [0 1]);colorbar
            end
    end
end
%%
% Optimizer
[Opt, Met]             = imregconfig('Multimodal');
Opt.Epsilon            = 1.5e-6;
Opt.MaximumIterations  = 300;
Opt.InitialRadius      = 6.25e-4;
registration_type = 'affine'; % 'translation', 'affine'

Mov = squeeze(B(:, :, 3));
Fxd = squeeze(A(:, :, 3));
RefRgt = imregtform(Mov, Fxd, registration_type, Opt, Met);
Mov_aft = imwarp(Mov,RefRgt,'OutputView',imref2d(size(Mov)));
%%
figure; imshow(Fxd(1:3000, 1:4000))
%%

figure;
subplot(1, 3, 1);
imshow(norm(double(Fxd)));
subplot(1, 3, 2);
imshow(norm(double(Mov)));
subplot(1, 3, 3);
imshow(norm(double(Mov_aft)));

%%
figure;
subplot(1, 2, 1)
% display_img = norm(double(Fxd - Mov_aft));
display_img = norm(double(Mov_aft)) - norm(double(Fxd));
colormap(parula)
imagesc(display_img); colorbar
box off
axis off
subplot(1, 2, 2)
display_img = zeros(size(Fxd, 1), size(Fxd, 2), 3);
display_img(:, :, 1) = norm(double(Fxd));
display_img(:, :, 2) = norm(double(Mov_aft));
% display_img(:, :, 3) = norm(double(Mov_aft));
imshow(display_img)