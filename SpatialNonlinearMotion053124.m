function OUT = SpatialNonlinearMotion053124(IN, OLED, DIO)

sw = OLED.width;
sh = OLED.height;

a2d = @(x) x*pi/180;
if IN.nDirection > 1
    IsDirecitonal = 1;
    Angs = round((0:IN.nDirection-1)*(360/IN.nDirection));
    DircVectors = [cos(a2d(Angs(:))), sin(a2d(Angs(:)))];
else
    IsDirecitonal = 0;
    DircVectors = [1, 0];
end

%% Calculate the size of the canvas
maxT = sum(IN.CycleTime(3:end, 1));
maxSpd = max(IN.Speeds)/OLED.pixelSize;
if IsDirecitonal
    maxW = ceil(maxT*maxSpd*2+sw);
    maxH = ceil(maxT*maxSpd*2+sh);
else
    maxW = ceil(maxT*maxSpd*2+sw);
    maxH = sh;
end
%%
pixel2um = OLED.pixelSize;
nStep = size(IN.CycleTime, 1);

%%
calibration_params = load('NF_IntensityMeasurement_032921_NF_2_002.mat', 'parameters');
calibration_params = calibration_params.parameters;

contrast_mean = ProjIntCalib(calibration_params, 0.5);
%% Canvas
SpatialWavelength = IN.TexureSpatialL/OLED.pixelSize;
nSpaPattern  = length(IN.TexureSpatialL) * IN.nTexure;
SpaPattern = nan(nSpaPattern, 5);
is_size_fixed = 1;
ii = 1;
for i = 1:length(IN.TexureSpatialL)
    for j = 1:IN.nTexure
        rand_seed_id = 20+j;
        SpaPattern(ii, :) = [i, j, SpatialWavelength(i), rand_seed_id, IN.TexureSpatialL(i)];
        c_canvas = gaussianblubs(maxW, maxH, SpatialWavelength(i), pixel2um, rand_seed_id, is_size_fixed);
        if ii == 1
            canvaset = nan(size(c_canvas, 2), size(c_canvas, 1), nSpaPattern);
        end
        canvaset(:, :, ii) = ProjIntCalib(calibration_params, c_canvas');
        ii = ii +1;
    end
end

canvasW = size(canvaset, 1);
canvasH = size(canvaset, 2);

%%  Mask
[X, Y] = meshgrid(1:size(canvaset, 2), 1:size(canvaset, 1));
X = X-0.5*(1+size(canvaset, 2));
Y = Y-0.5*(1+size(canvaset, 1));

bg_mask = sqrt(X.^2 + Y.^2) < IN.RFCenterSize*0.5/OLED.pixelSize;
% bg_mask = bg_mask';
% figure; imagesc(bg_mask); colorbar
%%
[SpaLs, Texures, Directs] = meshgrid(1:length(IN.TexureSpatialL), 1:IN.nTexure, 1:IN.nDirection);
Patterns = [SpaLs(:) Texures(:) Directs(:) DircVectors(Directs(:), :)];
nPattern = size(Patterns, 1);

%%
maskdia_adj = IN.MaskSize*0.5/OLED.pixelSize;
contrast_mean = 0.5;
%%
FrmTable = [];
FrmTime = [];
cgpencol(contrast_mean, contrast_mean, contrast_mean)
cgellipse(0,0,maskdia_adj,maskdia_adj,'f')
if ~OLED.IsDebug
    pin16high
end
OUT.startT = cgflip;

for i = 1:IN.nRepeat
    rng(42)
    PatternIds = randperm(nPattern);
    for k = 1:nPattern
        pid = PatternIds(k);

        vspd = IN.Speeds/OLED.pixelSize * Patterns(pid, 4:5);
        canvas = squeeze(canvaset(:, :, pid));
        cgloadarray(1, canvasW, canvasH, repmat(canvas(:), 1, 3), canvasW, canvasH);
        canvas_1 = canvas;
        canvas_1(bg_mask) = 1;
        canvas_2 = canvas;
        canvas_2(bg_mask) = 0;
        canvas_3 = canvas;
        canvas_3(bg_mask) = 0;

        bg_canvas = cat(3, canvas_1, canvas_2, canvas_3);

        cgloadarray(2, canvasW, canvasH, reshape(bg_canvas, [], 3), canvasW, canvasH);
        cgtrncol(2, 'r');
        clear canvas_1 canvas_2 canvas_3
        
        for q = 1:nStep
            cFrmTime = [];
            switch q
                case 1 % stationary
                    cgdrawsprite(1, 0, 0);
                    FlipTime = cgflip;
                    FrmTable = [FrmTable; FlipTime, i, k, q, pid];
                    cFrmTime = [cFrmTime; FlipTime, 0 0 nan(1, 4)];
                    tic
                    while toc < IN.CycleTime(q, 1) - 0.5*1/60
                    end
                    cgdrawsprite(1, 0, 0);
                    FlipTime = cgflip;
                    cFrmTime = [cFrmTime; FlipTime, 0 0 nan(1, 4)];
                case 2 % local motion
                    if IN.CycleTime(q, 1) == 0
                        pos = zeros(1, 2);
                        continue
                    end
                    tic
                    IsFirst = 1;
                    while toc < IN.CycleTime(q, 1) - 0.5*1/60
                        pos = toc*vspd;
                        cgdrawsprite(1, pos(1), pos(2));
                        cgdrawsprite(2, 0, 0);
                        FlipTime = cgflip;
                        if IsFirst
                            FrmTable = [FrmTable; FlipTime, i, k, q, pid];
                            IsFirst = 0;
                        end
                        cFrmTime = [cFrmTime; FlipTime, pos 0 0 nan(1, 2)];
                    end
                case 3 % surround motion
                    if IN.CycleTime(q, 1) == 0
                        pos = zeros(1, 2);
                        continue
                    end
                    tic
                    IsFirst = 1;
                    while toc < IN.CycleTime(q, 1) - 0.5*1/60
                        cpos = toc*vspd;
                        cgdrawsprite(1, cpos(1), cpos(2));
                        cgdrawsprite(3, 0, 0);
                        FlipTime = cgflip;
                        if IsFirst
                            FrmTable = [FrmTable; FlipTime, i, k, q, pid];
                            IsFirst = 0;
                        end
                        cFrmTime = [cFrmTime; FlipTime, cpos nan(1, 2) pos];
                    end
                case 4 % global motion
                    tic
                    IsFirst = 1;
                    while toc < IN.CycleTime(q, 1) - 0.5*1/60
                        cpos = pos+toc*vspd;
                        cgdrawsprite(1, cpos(1), cpos(2));

                        if IsFirst
                            FlipTime = cgflip;
                            FrmTable = [FrmTable; FlipTime, i, k, q, pid];
                            IsFirst = 0;
                        else
                            FlipTime = cgflip;
                        end
                        cFrmTime = [cFrmTime; FlipTime, cpos nan(1, 4)];
                    end
            end
            tic
            FrmTable = [FrmTable; FlipTime, i, k, q, pid];
            FrmTime = [FrmTime; cFrmTime];
            IsFirst = 1;
            while toc < IN.CycleTime(q, 2) - 0.5*1/60
                % general for surround motion
                if IsFirst && q == 2 && IN.CycleTime(3, 1) > 0
                    mbg_mask = sqrt((X+round(pos(2))).^2 + (Y+round(pos(1))).^2) < IN.RFCenterSize*0.5/OLED.pixelSize;
                    % mbg_mask = mbg_mask';
                    canvas_1 = ones(size(canvas));
                    canvas_1(bg_mask) = canvas(mbg_mask);
                    canvas_2 = zeros(size(canvas));
                    canvas_2(bg_mask) = canvas(mbg_mask);
                    canvas_3 = zeros(size(canvas));
                    canvas_3(bg_mask) = canvas(mbg_mask);

                    mbg_canvas = cat(3, canvas_1, canvas_2, canvas_3);
                    
                    cgloadarray(3, canvasW, canvasH, reshape(mbg_canvas, [], 3), canvasW, canvasH);
                    cgtrncol(3, 'r');
                    IsFirst = 0;
                    clear canvas_1 canvas_2 canvas_3
                end
            end
        end
    end
end
cgpencol(contrast_mean, contrast_mean, contrast_mean)
cgellipse(0,0,maskdia_adj,maskdia_adj,'f')
if ~OLED.IsDebug
    pin16low
end
OUT.endT = cgflip;
cgflip(0,0,0)
cgflip
OUT.FrmTable = FrmTable;
OUT.FrmTime = FrmTime;
OUT.PatternTable = Patterns;
OUT.SpaPattern = SpaPattern;
OUT.canvaset = canvaset; % with contrast correction