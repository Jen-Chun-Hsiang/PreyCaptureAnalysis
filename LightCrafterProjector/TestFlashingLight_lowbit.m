client = stage.core.network.StageClient();
client.connect();
canvas_size = client.getCanvasSize();
%%
target_frame_up_factor = 1; % 1, 2, 4, 6, 8, 12, 24
flashlight_propsize = 0.8; % from 0 to 1 (full screen)
stim_duration = 10;
max_led_current = 50;
is_antialias = 0;
is_prerender = 0;

%% Initiate LightCrafter modulation
addpath(genpath('C:\Users\kerschensteinerlab\Documents\MATLAB\matlab-lcr'))
refreshRate = client.getMonitorRefreshRate();
lightCrafter = LightCrafter4500(refreshRate);
lightCrafter.connect();
lightCrafter.setMode('pattern');
lightCrafter.setLedEnables(true, false, false, false);

%% Set maximum current
[led_currents(1), led_currents(2), led_currents(3)] = lightCrafter.getLedCurrents();
% [red_current, green_current, blue_current] = obj.lightCrafter.getLedCurrents();
lightCrafter.setLedCurrents(...
    min(led_currents(1), max_led_current),...
    min(led_currents(2), max_led_current),...
    min(led_currents(3), max_led_current));

%% 
patternRatesToAttributes = containers.Map('KeyType', 'double', 'ValueType', 'any');
patternRatesToAttributes(1 * refreshRate)  = {8, 'white', 1};
patternRatesToAttributes(2 * refreshRate)  = {8, 'white', 2};
patternRatesToAttributes(4 * refreshRate)  = {6, 'white', 4};
patternRatesToAttributes(6 * refreshRate)  = {4, 'white', 6};
patternRatesToAttributes(8 * refreshRate)  = {3, 'white', 8};
patternRatesToAttributes(12 * refreshRate) = {2, 'white', 12};
patternRatesToAttributes(24 * refreshRate) = {1, 'white', 24};

attributes = patternRatesToAttributes(target_frame_up_factor*refreshRate);
lightCrafter.setPatternAttributes(attributes{:});
renderer = stage.builtin.renderers.PatternRenderer(attributes{3}, attributes{1});
client.setCanvasRenderer(renderer);

current_patternrate = lightCrafter.currentPatternRate();

%%
flashlight = stage.builtin.stimuli.Rectangle();
flashlight.size = canvas_size*0.8;

%% flip between black and white
trackerColor = stage.builtin.controllers.PropertyController(flashlight, 'color', @(s)mod(s.frame, 2) && double(s.time + (1/s.frameRate) < presentation.duration));
presentation = stage.core.Presentation(stim_duration);
presentation.addStimulus(flashlight);
presentation.addController(trackerColor);  
% RENDER
if is_prerender
    player = stage.builtin.players.PrerenderedPlayer(presentation);
else
    player = stage.builtin.players.RealtimePlayer(presentation);
end
if ~is_antialias
    player.setCompositor(stage.builtin.compositors.PatternCompositor());
else
    player.setCompositor(sa_labs.util.ExactPatternCompositor());
end
client.play(player);