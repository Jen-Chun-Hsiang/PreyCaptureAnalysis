windowSize = [600, 800];

% Create the window in non-fullscreen mode
window = stage.core.Window(windowSize, false);
%window.open();

% Set the background to black
%window.setClearColor([0, 0, 0]);
% Create the canvas
canvas = stage.core.Canvas(window, 'disableDwm', false);

stageClient = stage.core.network.StageClient();
stageClient.connect('localhost', 5678);

numPatterns = 8;
bitDepth = 2;
renderer = stage.builtin.renderers.PatternRenderer(numPatterns, bitDepth);
stageClient.setCanvasRenderer(renderer);

presentationDuration = 2;
presentation = stage.core.Presentation(presentationDuration);

% Create a circle (spot) stimuls
circleRadius = 100;
circlePosition = windowSize/2;
cicle = stage.bultin.stimuli.Ellipse();
circle.color = [1, 1, 1];
circle.radiusX = circleRadius; 
circle.radiusY = circleRadius; 
circle.position = circlePosition; 

% Add the circle to the presentation
presentation.addStimulus(circle)

% Play the presentation
window.play(presentation);

window.close();
stageClient.disconnect();