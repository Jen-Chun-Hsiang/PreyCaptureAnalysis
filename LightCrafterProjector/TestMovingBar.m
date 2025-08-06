client = stage.core.network.StageClient();
client.connect();
canvas_size = client.getCanvasSize();
% Create the bar stimulus.
bar = stage.builtin.stimuli.Rectangle();
bar.size = [100, 100];
bar.color = [0, 0, 0];

% Create a controller to change the bar's position property as a function of time.
barPositionController = stage.builtin.controllers.PropertyController(bar, 'position', @(state)[state.time*110+100, canvas_size(1)/2]);

% Create a 5 second presentation and add the stimulus and controller.
presentation = stage.core.Presentation(5);
presentation.setBackgroundColor(1);
presentation.addStimulus(bar);
presentation.addController(barPositionController);

%p = stage.core.Presentation(5);
%p.addStimulus(stage.builtin.stimuli.Rectangle());
player = stage.builtin.players.RealtimePlayer(presentation);
client.play(player);