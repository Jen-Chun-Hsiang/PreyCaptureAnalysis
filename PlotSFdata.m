x = values(Ids{5});
y = values(Ids{6});
% p = ranksum(x,y)
[~, p] = ttest2(x, y) 

%%
figure; hold on
colormap(gray);
c = (stdSTA-min(stdSTA(:)))/range(stdSTA(:));
imagesc(400:600, 200:300, c');
plot(100+[0 100]/4.275, 200*ones(1, 2), 'w', 'LineWidth', 2);

%%
%--- (2) Show the image -------------------
I = stdSTA;
Imed = median( I(:) );
Imax = max(    I(:) );
a = 0.5 / (Imax - Imed);
b = 1   - a*Imax;
I = a*I + b;
I(I < 0) = 0;
I(I > 1) = 1;

I = I';
hFig = figure;
hAx  = axes(hFig);
imshow(I, 'Parent', hAx);
axis(hAx, "image");   % ensure 1:1 pixel aspect
hold(hAx, "on");

% (2) Define your pixel‐to‐unit conversion
%     e.g. if 1 pixel = 0.1 µm, then:
um_per_pixel = 4.257;

% (3) Decide the physical length of the scale bar
desiredBarLength_um = 100;              % want a 10-µm bar
barLength_px = round(desiredBarLength_um / um_per_pixel);

% (4) Figure out where to place the bar (in pixels)
[H, W, ~] = size(I);   % image height, width
margin = 10;           % leave 10 px from bottom & right edges

% Coordinates for the line: from (x_start, y_pos) to (x_end, y_pos)
x_end   = W - margin;                    % right end of bar
x_start = x_end - barLength_px;          % left end of bar
y_pos   = H - margin;                    % vertical position (just above bottom)

% (5) Draw the scale bar as a solid line
line([x_start, x_end], [y_pos, y_pos], ...
     'Color',     'white', ...
     'LineWidth', 2, ...
     'Parent',    hAx);

% (6) Add a label just above the bar, centered
x_text = (x_start + x_end) / 2;
y_text = y_pos - 5;   % 5 px above the line
text(x_text, y_text, sprintf("%d µm", desiredBarLength_um), ...
     'Color',               'white', ...
     'FontSize',            12, ...
     'HorizontalAlignment', 'center', ...
     'Parent',              hAx);


%%
% (1) Choose x (percent of each dimension)
x = 50;   % e.g. 50% → half width and half height

% (2) Get original size
[H, W, ~] = size(I);

% (3) Compute size of the cropped region
cropH = round(H * (x/100));
cropW = round(W * (x/100));

% (4) Compute start indices so that crop is centered
rowStart = floor((H - cropH)/2) + 1;
colStart = floor((W - cropW)/2) + 1;

% (5) Crop out the central region
Icrop = I( rowStart : rowStart + cropH - 1, ...
           colStart : colStart + cropW - 1, : );
hFig = figure;
hAx  = axes(hFig);
imshow(Icrop, 'Parent', hAx);
axis(hAx, "image");   % ensure 1:1 pixel aspect
hold(hAx, "on");

% (2) Define your pixel‐to‐unit conversion
%     e.g. if 1 pixel = 0.1 µm, then:
um_per_pixel = 4.257;

% (3) Decide the physical length of the scale bar
desiredBarLength_um = 100;              % want a 10-µm bar
barLength_px = round(desiredBarLength_um / um_per_pixel);

% (4) Figure out where to place the bar (in pixels)
[H, W, ~] = size(Icrop);   % image height, width
margin = 10;           % leave 10 px from bottom & right edges

% Coordinates for the line: from (x_start, y_pos) to (x_end, y_pos)
x_end   = W - margin;                    % right end of bar
x_start = x_end - barLength_px;          % left end of bar
y_pos   = H - margin;                    % vertical position (just above bottom)

% (5) Draw the scale bar as a solid line
line([x_start, x_end]-300, [y_pos, y_pos], ...
     'Color',     'white', ...
     'LineWidth', 2, ...
     'Parent',    hAx);

% (6) Add a label just above the bar, centered
x_text = (x_start + x_end) / 2-300;
y_text = y_pos - 15;   % 5 px above the line
text(x_text, y_text, sprintf("%d µm", desiredBarLength_um), ...
     'Color',               'white', ...
     'FontSize',            12, ...
     'HorizontalAlignment', 'center', ...
     'Parent',              hAx);



