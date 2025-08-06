classdef MovingBar < handle
    properties
        preTime = 250;
        stimTime = 5000; % Duration of the bar moving across the screen
        tailTime = 250;
        intensity = 1.0;
        barLength = 600; % microns
        barWidth = 200; % microns
        barSpeed = 1000; % microns per second
        distance = 3000; % microns
        angleOffset = 0;
        numberOfAngles = 12;
        numberOfCycles = 3;
    end
    
    methods
        function obj = MovingBar()
            % Constructor
        end
        
        function run(obj, lightCrafterDevice)
            % Setup presentation window
            window = stage.core.Window('fullScreen', false, 'position', [100, 100, 1024, 768]);
            angles = mod(round(0:360/obj.numberOfAngles:(360-.01)) + obj.angleOffset, 360);
            
            for cycle = 1:obj.numberOfCycles
                for angleIndex = 1:length(angles)
                    barAngle = angles(angleIndex);
                    presentation = obj.createPresentation(barAngle);
                    player = stage.core.Player(presentation, window);
                    player.play();
                end
            end
            
            window.close();
        end
        
        function presentation = createPresentation(obj, barAngle)
            presentationDuration = (obj.preTime + obj.stimTime + obj.tailTime) / 1000; % Convert ms to s
            presentation = stage.core.Presentation(presentationDuration);
            
            bar = stage.builtin.stimuli.Rectangle();
            bar.color = obj.intensity;
            bar.size = [obj.barLength, obj.barWidth]; % Conversion to pixels might be needed
            bar.orientation = barAngle;
            bar.position = [512, 384]; % Center of a 1024x768 window
            
            presentation.addStimulus(bar);
            
            % Movement based on speed and direction
            xStep = obj.barSpeed * cosd(barAngle) / 1000; % pixels per millisecond
            yStep = obj.barSpeed * sind(barAngle) / 1000; % pixels per millisecond
            
            movementController = stage.builtin.controllers.PropertyController(bar, 'position', ...
                @(state)[512 + xStep * (state.time - obj.preTime / 1000), ...
                         384 + yStep * (state.time - obj.preTime / 1000)]);
            presentation.addController(movementController);
        end
    end
end
