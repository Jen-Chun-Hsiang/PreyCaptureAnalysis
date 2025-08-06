client = stage.core.network.StageClient();
client.connect();
client.setMonitorGamma(1);

trueCanvasSize = client.getCanvasSize();
canvasSize = [trueCanvasSize(1) * 2, trueCanvasSize(2)];
client.setCanvasProjectionIdentity();
client.setCanvasProjectionOrthographic(0, canvasSize(1), 0, canvasSize(2));

%%
target_frame_up_factor = 2; % 1, 2, 4, 6, 8, 12, 24
flashlight_propsize = 0.8; % from 0 to 1 (full screen)
stim_duration = 10;
max_led_current = 100;
is_antialias = 0;
is_prerender = 0;
LED_color_str = 'white'; % 'none', 'red', 'green', 'yellow'
%% Initiate LightCrafter modulation
addpath(genpath('C:\Users\kerschensteinerlab\Documents\MATLAB\matlab-lcr'))
refreshRate = client.getMonitorRefreshRate();
lightCrafter = LightCrafter4500(refreshRate);
lightCrafter.connect();
lightCrafter.setMode('pattern');
lightCrafter.setLedEnables(true, false, false, false);
 [auto, red, green, blue] = lightCrafter.getLedEnables();

%% Set maximum current
[led_currents(1), led_currents(2), led_currents(3)] = lightCrafter.getLedCurrents();
% [red_current, green_current, blue_current] = obj.lightCrafter.getLedCurrents();
lightCrafter.setLedCurrents(...
    min(led_currents(1), max_led_current),...
    min(led_currents(2), max_led_current),...
    min(led_currents(3), max_led_current));

%% 
patternRatesToAttributes = containers.Map('KeyType', 'double', 'ValueType', 'any');

patternRatesToAttributes(1 * refreshRate)  = {8, LED_color_str, 1};
patternRatesToAttributes(2 * refreshRate)  = {8, LED_color_str, 2};
patternRatesToAttributes(4 * refreshRate)  = {6, LED_color_str, 4};
patternRatesToAttributes(6 * refreshRate)  = {4, LED_color_str, 6};
patternRatesToAttributes(8 * refreshRate)  = {3, LED_color_str, 8};
patternRatesToAttributes(12 * refreshRate) = {2, LED_color_str, 12};
patternRatesToAttributes(24 * refreshRate) = {1, LED_color_str, 24};

attributes = patternRatesToAttributes(target_frame_up_factor*refreshRate);
lightCrafter.setPatternAttributes(attributes{:});
            
renderer = stage.builtin.renderers.PatternRenderer(attributes{3}, attributes{1});
client.setCanvasRenderer(renderer);

current_patternrate = lightCrafter.currentPatternRate();

%% Create the bar stimulus.
bar = stage.builtin.stimuli.Rectangle();
bar.size = [100, 100];
bar.color = 1*ones(1, 3);
bar.position = canvasSize*0.5;

% Create a controller to change the bar's position property as a function of time.
% barPositionController = stage.builtin.controllers.PropertyController(bar,...
%     'position', @(state)[state.time*100+100, canvasSize(1)*0.5]);

barPositionController = stage.builtin.controllers.PropertyController(bar,...
    'position', @(state)[canvasSize(2), state.time*1000+100]);

% [X, Y] = meshgrid(-600:200:600, -600:200:600);
% all_position = [X(:), Y(:)];
% num_position = size(all_position, 1);
% all_position = all_position(randperm(num_position), :);
% barPositionController = stage.builtin.controllers.PropertyController(bar, 'position',...
%     @(state)[all_position(ceil(mod(state.time, num_position)+1), 1),...
%              all_position(ceil(mod(state.time, num_position)+1), 2)]);
%%
presentation = stage.core.Presentation(stim_duration);
presentation.setBackgroundColor(0);
presentation.addStimulus(bar);
presentation.addController(barPositionController);  
% RENDER
if is_prerender
    player = stage.builtin.players.PrerenderedPlayer(presentation);
else
    player = stage.builtin.players.RealtimePlayer(presentation);
end
player.setCompositor(stage.builtin.compositors.PatternCompositor());

% if ~is_antialias
%     player.setCompositor(stage.builtin.compositors.PatternCompositor());
% else
%     player.setCompositor(sa_labs.util.ExactPatternCompositor());
% end
client.play(player);
play_info = client.getPlayInfo();
disp('presentation ends')