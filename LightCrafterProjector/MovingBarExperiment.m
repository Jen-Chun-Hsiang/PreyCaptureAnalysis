classdef MovingBarExperiment < BaseProtocol
    methods
        function obj = MovingBarExperiment()
            % Constructor
            obj@BaseProtocol(); % Call BaseProtocol constructor
        end
        
        function run(obj)
            % Specify settings for LightCrafter
            settings = struct();
            settings.refreshRate = 120; % Example: 120 Hz
            settings.colorMode = 'standard'; % or 'uv' based on your setup
            settings.bitDepth = 8; % Example: 8-bit depth for higher frame rate
            
            lightCrafterDevice = LightCrafterDevice(settings);
            movingBar = MovingBar();
            
             obj.prepareRun();
            for i = 1:obj.totalNumEpochs
                disp(['Epoch ' num2str(i)]);
                movingBar.run(lightCrafterDevice); % Display the moving bar
            end
            obj.completeRun();
            
            lightCrafterDevice.close(); % Ensure device is properly disconnected
        end
    end
end