classdef DriftingGratings < handle
    properties
        preTime = 250; % ms
        tailTime = 250; % ms
        stimTime = 5000; % ms
        
        movementDelay = 200; % ms
        
        gratingWidth = 3000; % um
        gratingLength = 3000; % um
        gratingSpeed = 1200; % um/s
        cycleHalfWidth = 300; % um
        apertureDiameter = 0; % um
        gratingProfile = 'square'; % 'sine', 'square', or 'sawtooth'
        contrast = 1;
        
        numberOfAngles = 12;
        numberOfCycles = 2;
    end
    
    properties (Dependent)
        spatialFreq % cycles/degree
        temporalFreq % cycles/s (Hz)
        totalNumEpochs
    end
    
    properties (Access = private)
        curAngle
        angles
    end
    
    methods
        function obj = DriftingGratings()
            % Initialize angles
            obj.angles = mod(0:360/obj.numberOfAngles:359, 360);
        end
        
        function prepareRun(obj)
            % Randomize angles for each run
            obj.angles = obj.angles(randperm(obj.numberOfAngles));
        end
        
        function [p, grat] = createPresentation(obj, angle)
            % Create the drifting grating presentation
            obj.curAngle = angle;
            
            % Assuming Psychtoolbox is used for presentation
            p = []; % Placeholder for presentation object
            grat = []; % Placeholder for grating object
            
            % Implementation details for creating and configuring the
            % drifting grating stimulus would go here. This could involve
            % configuring a Psychtoolbox texture and drawing it at various
            % orientations and phases to achieve the drifting effect.
        end
        
        function run(obj)
            % Main method to execute the protocol
            for cycle = 1:obj.numberOfCycles
                for angleIndex = 1:length(obj.angles)
                    [p, grat] = obj.createPresentation(obj.angles(angleIndex));
                    % Present the drifting grating here
                    % This would involve rendering the presentation `p`
                    % created by `createPresentation`
                end
            end
        end
        
        function totalNumEpochs = get.totalNumEpochs(obj)
            totalNumEpochs = obj.numberOfCycles * obj.numberOfAngles;
        end
        
        function spatialFreq = get.spatialFreq(obj)
            % Calculate spatial frequency based on cycleHalfWidth
            micronperdeg = 30; % Assuming 1 deg visual angle = 30um (mouse retina)
            spatialFreq = 1 / (obj.cycleHalfWidth * 2 / micronperdeg);
        end
        
        function temporalFreq = get.temporalFreq(obj)
            % Calculate temporal frequency based on gratingSpeed and cycleHalfWidth
            temporalFreq = obj.gratingSpeed / (obj.cycleHalfWidth * 2);
        end
    end
end
