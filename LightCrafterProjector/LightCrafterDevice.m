classdef LightCrafterDevice < handle
    properties (Access = private)
        stageClient
        lightCrafter  % Instance of LightCrafter4500
%         orientation
%         baseTranslation = [0,0]
    end
    
    methods
        function obj = LightCrafterDevice(settings)
            %% Initialize with provided settings or defaults
            if nargin < 1 || isempty(settings)
                settings = obj.defaultSettings();
            end
            
            %% Setup Stage and LightCrafter connection
            obj.stageClient = stage.core.network.StageClient();
            obj.stageClient.connect(settings.host, settings.port);
            obj.stageClient.setMonitorGamma(1);
            
            canvasSize = obj.calculateCanvasSize(obj.stageClient.getCanvasSize());
            obj.stageClient.setCanvasProjectionIdentity();
            obj.stageClient.setCanvasProjectionOrthographic(0, canvasSize(1), 0, canvasSize(2));
            
            %% LightCrafter setup
%             obj.orientation = settings.orientation;
            % Initialize LightCrafter4500 with specified refresh rate and color mode
            obj.lightCrafter = LightCrafter4500(settings.refreshRate, settings.colorMode);
            
            % Connect to the LightCrafter
            obj.lightCrafter.connect();
            
            % Configure bit depth and refresh rate
            obj.configureDevice(settings);
%             obj.baseTranslation = settings.canvasTranslation;
            
            %% Apply settings
            obj.applySettings(settings);
        end
        
        function settings = defaultSettings()
            %% Default configuration settings
            settings = struct();
            settings.host = 'localhost';
            settings.port = 5678;
            settings.projectorColorMode = 'standard';
            settings.bitDepth = 8; % Lower bit depth for higher frame rate
            settings.frameRate = 120; % Higher frame rate
            
            settings.orientation = [0,0];
            settings.micronsPerPixel = 2;
            settings.canvasTranslation = [0,0];
            settings.angleOffset = 0;
            % Add more default settings as needed
        end
        
        function configureDevice(obj, settings)
            % Set mode to pattern display for lower bit depth and higher refresh rate
            obj.lightCrafter.setMode('pattern');
            
            % Configure LED enables, if required
            % obj.lightCrafter.setLedEnables(true, true, true, true); % Example: Enable all LEDs
            
            % Configure the pattern sequence attributes
            color = {'red'}; % Example: Use red LED for simplicity. Modify as needed.
            obj.lightCrafter.setPatternAttributes(settings.bitDepth, color, []);
        end
        
        function canvasSize = calculateCanvasSize(~, trueCanvasSize)
            %% Calculate adjusted canvas size if needed
            canvasSize = [trueCanvasSize(1) * 2, trueCanvasSize(2)];
        end
        
        function lcr = setupLightCrafter(obj, settings)
            %% Setup LightCrafter connection and configuration
            % Placeholder for creating and configuring LightCrafter instance
            % This method should establish a connection to the LightCrafter and configure it according to 'settings'
            disp('LightCrafter configured with settings:');
            disp(settings);
            lcr = []; % Replace with actual LightCrafter setup code
        end
        
        function applySettings(obj, settings)
            %% Apply various settings to LightCrafter and stage client
            % This method applies the settings such as orientation, color mode, etc., to the LightCrafter and stage projection
        end
        
        function close(obj)
            %% Cleanup before closing the object
            if ~isempty(obj.stageClient)
                obj.stageClient.disconnect();
            end
            if ~isempty(obj.lightCrafter)
                % Disconnect LightCrafter here
            end
        end
    end
end
