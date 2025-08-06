classdef LightCrafterControl < handle
    properties (Access = private)
        % Assume these are handles or identifiers to interface with the actual LightCrafter device.
        lightCrafter
        ledEnablesCheckboxes = struct()
        patternRatePopupMenu
    end

    methods
        function obj = LightCrafterControl()
            % Constructor - Create the UI and initialize the LightCrafter connection
            obj.createUi();
            obj.initializeLightCrafter();
            obj.populateLedEnablesCheckboxes();
            obj.populatePatternRateList();
        end
        
        function createUi(obj)
            % Create the MATLAB UI components for controlling the LightCrafter
            
            figureHandle = figure('Name', 'LightCrafter Control', ...
                                  'Position', [100, 100, 320, 75], ...
                                  'MenuBar', 'none', ...
                                  'NumberTitle', 'off', ...
                                  'Resize', 'off');
            mainLayout = uix.HBox('Parent', figureHandle, 'Padding', 11, 'Spacing', 7);

            lightCrafterLayout = uix.Grid('Parent', mainLayout, 'Spacing', 7);
            uicontrol('Parent', lightCrafterLayout, 'Style', 'text', 'String', 'LED enables:');
            uicontrol('Parent', lightCrafterLayout, 'Style', 'text', 'String', 'Pattern rate:');

            ledEnablesLayout = uix.HBox('Parent', lightCrafterLayout, 'Spacing', 3);
            fields = {'auto', 'red', 'green', 'blue'};
            for i = 1:length(fields)
                obj.ledEnablesCheckboxes.(fields{i}) = uicontrol('Parent', ledEnablesLayout, ...
                                                                'Style', 'checkbox', ...
                                                                'String', ucfirst(fields{i}), ...
                                                                'Callback', @(src, evt) obj.onSelectedLedEnable());
            end

            obj.patternRatePopupMenu = uicontrol('Parent', lightCrafterLayout, ...
                                                 'Style', 'popupmenu', ...
                                                 'String', {' '}, ...
                                                 'Callback', @(src, evt) obj.onSelectedPatternRate());

            set(lightCrafterLayout, 'Widths', [70 -1], 'Heights', [23 23]);
        end
    end

    methods (Access = private)
        function initializeLightCrafter(obj)
            % Initialize communication with the LightCrafter device
            % This is where you would establish a connection to the device
            % and assign it to obj.lightCrafter.
        end

        function populateLedEnablesCheckboxes(obj)
            % Query the LightCrafter for the current LED enable states and update the UI
            % Example static values, replace these with actual queries to the device
            ledStates = struct('auto', 0, 'red', 1, 'green', 1, 'blue', 1);
            fields = fieldnames(obj.ledEnablesCheckboxes);
            for i = 1:length(fields)
                set(obj.ledEnablesCheckboxes.(fields{i}), 'Value', ledStates.(fields{i}));
            end
        end

        function onSelectedLedEnable(obj)
            % Callback for LED enable checkboxes
            % Update the LightCrafter's LED enable state based on the UI
        end

        function populatePatternRateList(obj)
            % Populate the pattern rate popup menu with available rates
            % Example static rates, replace these with actual queries to the device
            rates = {'60 Hz', '120 Hz', '240 Hz'};
            set(obj.patternRatePopupMenu, 'String', rates);
        end

        function onSelectedPatternRate(obj)
            % Callback for the pattern rate popup menu
            % Update the LightCrafter's pattern rate based on the UI selection
        end
    end

    methods (Access = private, Static)
        function str = ucfirst(str)
            % Utility function to capitalize the first letter of a string
            str = regexprep(str, '^(.)', '${upper($1)}');
        end
    end
end
