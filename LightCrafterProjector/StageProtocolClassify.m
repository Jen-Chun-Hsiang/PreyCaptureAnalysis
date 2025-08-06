classdef (Abstract) StageProtocolClassify < handle
    properties
        meanLevel = 0.0               % Background light intensity (0-1)
        meanLevel1 = 0.5              % Background intensity value pattern 1
        meanLevel2 = 0.5              % Background intensity value pattern 2
        contrast1 = 1                 % Weber contrast from mean for object, color 1
        contrast2 = 1                 % Weber contrast from mean for object, color 2
        offsetX = 0                   % um
        offsetY = 0                   % um
        
        NDF = 5                       % Filter NDF value
        blueLED = 24                  % 0-255
        greenLED = 0                  % 0-255
        redLED = 0                    % 0-255
        uvLED = 0
        colorPattern1 = 'blue'
        colorPattern2 = 'none'
        colorPattern3 = 'none'
        primaryObjectPattern = 1
        secondaryObjectPattern = 1
        backgroundPattern = 2
        colorCombinationMode = 'contrast'
    end
    
    properties (Dependent)
        numberOfPatterns              % Number of color patterns
        bitDepth = 8                  % Color bit depth
        prerender = false             % Enable prerendering to reduce frame dropping
        frameRate = 60                % Frame rate, not implemented for change
    end
    
    methods (Abstract)
        p = createPresentation(obj);  % Abstract method for creating a presentation
    end
    
    methods
        function obj = StageProtocolClassify()
            % Constructor
        end
        
        function numberOfPatterns = get.numberOfPatterns(obj)
            if ~strcmp(obj.colorPattern3, 'none')
                numberOfPatterns = 3;
            elseif ~strcmp(obj.colorPattern2, 'none')
                numberOfPatterns = 2;
            else
                numberOfPatterns = 1;
            end
        end
        
        function bitDepth = get.bitDepth(obj)
            if obj.numberOfPatterns < 3
                bitDepth = 8;
            else
                bitDepth = 6; % Adjust based on actual implementation constraints
            end
        end
        
        function prerender = get.prerender(obj)
            % Determine if prerendering is necessary based on the number of patterns
            prerender = obj.numberOfPatterns > 1;
        end
        
        function frameRate = get.frameRate(obj)
            % Frame rate is fixed for this example
            frameRate = 60;
        end
    end
    
    methods (Access = protected)
        function [pround, p] = um2pix(obj, um)
            % Example conversion function from micrometers to pixels
            micronsPerPixel = 1; % This should be defined or calculated based on your setup
            p = um / micronsPerPixel;
            pround = round(p);
        end
    end
end
