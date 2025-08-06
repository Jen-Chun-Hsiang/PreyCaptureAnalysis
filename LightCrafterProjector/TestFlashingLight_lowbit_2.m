client = stage.core.network.StageClient();
client.connect();
client.setMonitorGamma(1);

trueCanvasSize = client.getCanvasSize();
canvasSize = [trueCanvasSize(1) * 2, trueCanvasSize(2)];
client.setCanvasProjectionIdentity();
client.setCanvasProjectionOrthographic(0, canvasSize(1), 0, canvasSize(2));

%%
target_frame_up_factor = 12; % 1, 2, 4, 6, 8, 12, 24
speed_fac = 30;
flashlight_propsize = 1; % from 0 to 1 (full screen)
stim_duration = 10;
max_led_current = 10;
is_antialias = 0;
is_prerender = 0;
LED_color_str = 'green'; % 'white', 'none', 'red', 'green', 'yellow'
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

%%
flashlight = stage.builtin.stimuli.Rectangle();
flashlight.size = canvasSize*1;
flashlight.position = canvasSize*0.5;



%% flip between black and white
% trackerColor = stage.builtin.controllers.PropertyController(flashlight, 'color',...
%     @(s)mod(s.frame, 2)>=1);
trackerColor = stage.builtin.controllers.PropertyController(flashlight, 'color',...
    @(s)mod(s.time*speed_fac, 2)>1);
% trackerColor = stage.builtin.controllers.PropertyController(flashlight, 'color',...
%     @(s)(mod(s.frame*speed_fac, 2)>1) && double(s.time + (1/s.frameRate) < presentation.duration));
% (mod(s.frame*speed_fac, 2)>1) 
%%
presentation = stage.core.Presentation(stim_duration);
%presentation.setBackgroundColor(1);
presentation.addStimulus(flashlight);
presentation.addController(trackerColor);  
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