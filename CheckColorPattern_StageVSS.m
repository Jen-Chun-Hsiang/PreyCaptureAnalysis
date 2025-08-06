intensity = 0:0.1:1;
colors = repmat(intensity', 1, 3);
nData = size(colors, 1);
Data = nan(nData, 4);
for i = 1:nData
    color = colors(i, :);
    numPattern =12;
    patternBitDepth = 2;
    patternIndex = 0;
    colorBitDepth = 8;

    c = mean(color(1:3));
    if c > 1
        c = 1;
    elseif c < 0
        c = 0;
    end
    patternColor = round(c * (2^patternBitDepth - 1));
    patternColor = bitshift(patternColor, patternIndex * patternBitDepth);

    % Split shifted pattern color into GRB components.
    bitMask = 2 ^ colorBitDepth - 1;
    g = bitand(bitshift(patternColor, -0 * colorBitDepth), bitMask);
    r = bitand(bitshift(patternColor, -1 * colorBitDepth), bitMask);
    b = bitand(bitshift(patternColor, -2 * colorBitDepth), bitMask);

    % Normalize and combine.
    g = g / bitMask;
    r = r / bitMask;
    b = b / bitMask;
    color = [r, g, b];
    Data(i, :) = [colors(i, 1) color];
end