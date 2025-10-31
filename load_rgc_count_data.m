% load_rgc_count_data.m
% Script to load RGC count data from Excel file
%
% Author: GitHub Copilot
% Date: October 30, 2025

% Clear workspace
clear;
clc;

% Define the folder path and file name
folderPath = '\\storage1.ris.wustl.edu\kerschensteinerd\Active\Emily\PreyCaptureRGC\rgc_count';
fileName = 'CellCounter_102725_11_11_headmodified.xlsx';

% Create full file path
filePath = fullfile(folderPath, fileName);

% Check if file exists
if ~isfile(filePath)
    error('File not found: %s', filePath);
end

% Load the Excel file into a table
try
    fprintf('Loading data from: %s\n', filePath);
    dataTable = readtable(filePath);
    
    % Display table information
    fprintf('\nTable loaded successfully!\n');
    fprintf('Number of rows: %d\n', height(dataTable));
    fprintf('Number of columns: %d\n', width(dataTable));
    fprintf('\nColumn names:\n');
    disp(dataTable.Properties.VariableNames);
    
    % Display first few rows
    fprintf('\nFirst few rows of data:\n');
    disp(head(dataTable));
    
    % Display data summary
    fprintf('\nData summary:\n');
    summary(dataTable);
    
catch ME
    error('Error loading Excel file: %s\nError message: %s', filePath, ME.message);
end

% The data is now available in the variable 'dataTable'
fprintf('\nData is available in the variable ''dataTable''\n');

%% Process the data table
fprintf('\n--- Processing data ---\n');

% (1) Remove the 'FileName' column
if ismember('FileName', dataTable.Properties.VariableNames)
    dataTable.FileName = [];
    fprintf('Removed ''FileName'' column\n');
end

% (2) Split data into two groups based on Marker_type
type1_data = dataTable(dataTable.Marker_type == 1, :);
type2_data = dataTable(dataTable.Marker_type == 2, :);

fprintf('Marker_type = 1: %d points\n', height(type1_data));
fprintf('Marker_type = 2: %d points\n', height(type2_data));

%% Prepare adaptive KDE heatmap
boundaryX = [];
boundaryY = [];
x_boundary = [];
y_boundary = [];
if height(type2_data) > 0
    boundaryX = type2_data.Marker_x;
    boundaryY = type2_data.Marker_y;
    x_boundary = [boundaryX; boundaryX(1)];
    y_boundary = [boundaryY; boundaryY(1)];
end

hasBoundary = height(type2_data) > 2;
hasCells = height(type1_data) > 0;
densityGrid = [];
xVec = [];
yVec = [];
hasDensity = false;

if hasCells && hasBoundary
    cellPoints = [type1_data.Marker_x, type1_data.Marker_y];
    try
        [densityGrid, xVec, yVec] = computeAdaptiveDensity(cellPoints, boundaryX, boundaryY);
        hasDensity = ~isempty(densityGrid);
        if hasDensity
            fprintf('Computed adaptive KDE density inside retina boundary.\n');
        else
            fprintf('Adaptive KDE skipped: insufficient evaluation points.\n');
        end
    catch kdeErr
        warning('AdaptiveKDE:Failure', 'Adaptive KDE failed: %s', kdeErr.message);
    end
end

%% Create overlay plot with KDE heatmap
figure;
hold on;

if hasDensity
    densityAlpha = 0.65 * double(~isnan(densityGrid));
    hImg = imagesc(xVec, yVec, densityGrid);
    set(hImg, 'AlphaData', densityAlpha);
    colormap(gca, hot);
    colorbar;
end

if hasCells
    scatter(type1_data.Marker_x, type1_data.Marker_y, 5, 'k', 'filled');
    fprintf('Plotted %d points (Marker_type = 1) as black dots\n', height(type1_data));
end

if height(type2_data) > 0
    plot(x_boundary, y_boundary, 'b-', 'LineWidth', 1);
    fprintf('Plotted boundary (Marker_type = 2) as blue line\n');
end

hold off;

set(gca, 'YDir', 'normal');
xlabel('X');
ylabel('Y');
title('RGC Count Data - Points, Boundary, Adaptive KDE');
axis equal;
grid on;

if hasDensity
    fprintf('Generated overlay plot with adaptive KDE heatmap.\n');
end

%% Standalone heatmap figure
if hasDensity
    figure;
    hHeat = imagesc(xVec, yVec, densityGrid);
    set(hHeat, 'AlphaData', double(~isnan(densityGrid)));
    set(gca, 'YDir', 'normal');
    colormap(gca, hot);
    colorbar;
    hold on;
    if height(type2_data) > 0
        plot(x_boundary, y_boundary, 'w-', 'LineWidth', 1);
    end
    hold off;
    axis equal;
    xlabel('X');
    ylabel('Y');
    title('Adaptive KDE Heatmap (Retina Boundary)');
    fprintf('Generated standalone adaptive KDE heatmap figure.\n');
end

fprintf('\nProcessing complete.\n');

function [densityGrid, xVec, yVec] = computeAdaptiveDensity(points, boundaryX, boundaryY)
% computeAdaptiveDensity - Evaluate adaptive KDE on a grid inside a polygon

    numPoints = size(points, 1);
    if numPoints < 5
        densityGrid = [];
        xVec = [];
        yVec = [];
        return;
    end

    minX = min(boundaryX);
    maxX = max(boundaryX);
    minY = min(boundaryY);
    maxY = max(boundaryY);

    gridResolution = max(120, min(250, round(sqrt(numPoints) * 12)));
    xVec = linspace(minX, maxX, gridResolution);
    yVec = linspace(minY, maxY, gridResolution);
    [gridX, gridY] = meshgrid(xVec, yVec);

    insideMask = inpolygon(gridX, gridY, boundaryX, boundaryY);
    evalPoints = [gridX(insideMask), gridY(insideMask)];

    if isempty(evalPoints)
        densityGrid = [];
        return;
    end

    densityValues = adaptiveKDE(points, evalPoints);

    densityGrid = nan(size(gridX));
    densityGrid(insideMask) = densityValues;
end

function densityVals = adaptiveKDE(points, evalPoints)
% adaptiveKDE - Adaptive Gaussian KDE with Abramson's method

    n = size(points, 1);
    if n < 2
        densityVals = zeros(size(evalPoints, 1), 1);
        return;
    end

    sigmaBase = std(points, 0, 1);
    sigmaBase(isnan(sigmaBase)) = 0;
    rangeVals = max(points, [], 1) - min(points, [], 1);
    zeroIdx = sigmaBase <= 0;
    sigmaBase(zeroIdx) = rangeVals(zeroIdx) / max(sqrt(n), 1);
    sigmaBase(sigmaBase <= 0) = max(rangeVals) * 0.01 + eps;

    baseBandwidth = 1.06 .* sigmaBase .* n^(-1/6);
    baseBandwidth(baseBandwidth <= 0) = eps;

    pilot = mvksdensity(points, points, 'Bandwidth', baseBandwidth, 'Kernel', 'normal');
    pilot = max(pilot, realmin);
    geomMeanPilot = exp(mean(log(pilot)));
    lambda = (pilot / geomMeanPilot) .^ (-0.5);

    densityVals = zeros(size(evalPoints, 1), 1);
    for i = 1:n
        hi = baseBandwidth .* lambda(i);
        hi(hi <= 0) = eps;
        dx = (evalPoints(:, 1) - points(i, 1)) ./ hi(1);
        dy = (evalPoints(:, 2) - points(i, 2)) ./ hi(2);
        kernelVals = exp(-0.5 .* (dx.^2 + dy.^2)) ./ (2 * pi * hi(1) * hi(2));
        densityVals = densityVals + kernelVals;
    end

    densityVals = densityVals / n;
end
