classdef (Abstract) StageProtocol
    properties
        meanLevel = 0.0     % Background light intensity (0-1)
        contrast = 1        % Weber contrast from mean for object
        offsetX = 0         % um
        offsetY = 0         % um
        NDF = 3             % Filter NDF value
        colorMode = 'standard' % 'standard', 'uv', etc.
        imaging = false     % Declare whether this is an imaging trial
        imagingFieldWidth = 125
        imagingFieldHeight = 125
        % Add more properties as needed
    end

    properties (Dependent)
        % Define dependent properties that calculate values based on other properties
    end

    methods (Abstract)
        p = createPresentation(obj); % Abstract method to create a visual presentation
    end
    
    methods
        function obj = StageProtocol()
            % Constructor for the protocol
            % Initialize any properties or state here
        end

        function run(obj)
            % Method to execute the protocol
            window = obj.setupWindow(); % Setup the presentation window
            presentation = obj.createPresentation(); % Create the presentation
            player = stage.core.Player(presentation, window);
            player.play(); % Play the presentation
            window.close(); % Close the window when done
        end

        function window = setupWindow(obj)
            % Setup the presentation window
            window = stage.core.Window('fullScreen', false, 'position', [100, 100, 1024, 768]);
            % Adjust window properties as needed
        end
    end

    methods (Access = protected)
        function [pround, p] = um2pix(obj, um)
            % Example conversion function from micrometers to pixels
            micronsPerPixel = 1; % Define the conversion rate based on your setup
            p = um / micronsPerPixel;
            pround = round(p);
        end
    end
end
