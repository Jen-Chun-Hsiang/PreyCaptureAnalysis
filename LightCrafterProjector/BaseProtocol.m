classdef (Abstract) BaseProtocol < handle
    properties
        chan1 = 'Amp1';
        chan1Mode = 'Cell attached';
        chan1Hold = 0;
        
        chan2 = 'None';
        chan2Mode = 'Off';
        chan2Hold = 0;
        
        chan3 = 'None';
        chan3Mode = 'Off';
        chan3Hold = 0;
        
        chan4 = 'None';
        chan4Mode = 'Off';
        chan4Hold = 0;
        
        spikeDetectorMode = 'advanced';
        spikeThreshold = -6; % pA or (pseudo-)std
        
        scanHeadTrigger = false; % Scanhead trigger for functional imaging
        stimTimeRecord = true; % Record stimulation time
        
        % Abstract properties to be defined by subclasses
        preTime;
        stimTime;
        tailTime;
        responsePlotMode;
        totalNumEpochs;
    end
    
    methods (Abstract)
        % Abstract method to be implemented by subclasses for creating the presentation
        preparePresentation(obj);
    end
    
    methods
        function obj = BaseProtocol()
            % Constructor for the protocol class
        end
        
        function prepareRun(obj)
            % Prepare the protocol run
            % This method can be overridden or extended by subclasses
            disp('Preparing protocol run...');
            % Initialize or reset states as needed
        end
        
        function prepareEpoch(obj)
            % Prepare individual epochs for the protocol
            disp('Preparing epoch...');
            % Setup for a single epoch; configure stimuli, responses, etc.
        end
        
        function run(obj)
            % Main method to execute the protocol
            obj.prepareRun();
            for i = 1:obj.totalNumEpochs
                disp(['Starting epoch ' num2str(i) ' of ' num2str(obj.totalNumEpochs)]);
                obj.prepareEpoch();
                % Simulate or perform epoch-specific actions here
                % For example, collecting data, showing stimuli, etc.
                pause(obj.stimTime / 1000); % Simulate the stimulation time
                disp(['Completed epoch ' num2str(i)]);
            end
            obj.completeRun();
        end
        
        function completeRun(obj)
            % Cleanup after the protocol run is complete
            disp('Protocol run completed.');
        end
        
        % Other methods related to device setup, epoch preparation, and cleanup can be added here
    end
end
