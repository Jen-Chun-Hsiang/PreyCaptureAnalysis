close all

client = stage.core.network.StageClient();
client.connect();
canvas_size = client.getCanvasSize();

% Create the flash stimulus.
flashlight = stage.builtin.stimuli.Rectangle();
flashlight.size = canvas_size*0.8;
%flashlight.color = [0, 0, 0];

flashlight.position = canvas_size*0.5;


% Create a controller to change the bar's position property as a function of time.
%barPositionController = stage.builtin.controllers.PropertyController(flashlight, 'position', @(state)[state.time*110+100, canvas_size(1)/2]);

% trackerColor = stage.builtin.controllers.PropertyController(flashlight, 'color',...
%     @(s)(mod(s.frame*1, 2)<1) && double(s.time + (1/s.frameRate) < presentation.duration));
trackerColor = stage.builtin.controllers.PropertyController(flashlight, 'color',...
    @(s)mod(s.time*30, 2)>1);

% Create a 5 second presentation and add the stimulus and controller.
presentation = stage.core.Presentation(5);
%presentation.setBackgroundColor(1);
presentation.addStimulus(flashlight);
presentation.addController(trackerColor);  
%presentation.addController(barPositionController);

%p = stage.core.Presentation(5);
%p.addStimulus(stage.builtin.stimuli.Rectangle());
% player = stage.builtin.players.RealtimePlayer(presentation);
% player = stage.builtin.players.PrerenderedPlayer(presentation);
time_thread = tic;
player = stage.builtin.players.TimeRecordedPrerenderPlayer(presentation, time_thread);
client.play(player);
play_info_1 = client.getPlayInfo();
player = stage.builtin.players.PrerenderedPlayer(presentation);
tic
client.play(player);
play_info = client.getPlayInfo();

%%
x = play_info.flipDurations;
y = diff(play_info_1.flip_times);

%%
figure; 
subplot(1, 2, 1); hold on
scatter(x, y, 5, 'k', 'filled');
xlabel('(built in) frame flip time (s)')
ylabel('(tic toc) frame flip time (s)');
plot([0.015 0.02], [0.015 0.02], '--b');
subplot(1, 2, 2); 
plot(y(:)./x(:), 'k');
xlabel('Epoch');
ylabel('(tic toc) / (built in)');
box off