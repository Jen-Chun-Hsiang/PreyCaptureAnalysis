function varargout  = WNAlign_sliderGUI(IN,OLED)
    bar_length = 200;
    % Create the main figure window
    fig = uifigure('Name', 'Slider GUI', 'Position', [100 100 550 400]);

    % Create sliders for Channel 1
    hSlider1 = uislider(fig, 'Orientation', 'horizontal', 'Position', [50 100 bar_length 3], 'Limits', [-100 100]);
    vSlider1 = uislider(fig, 'Orientation', 'vertical',   'Position', [50 120 3 bar_length], 'Limits', [-100 100]);

    % Create sliders for Channel 2
    hSlider2 = uislider(fig, 'Orientation', 'horizontal', 'Position', [300 100 bar_length 3], 'Limits', [-100 100]);
    vSlider2 = uislider(fig, 'Orientation', 'vertical',   'Position', [300 120 3 bar_length], 'Limits', [-100 100]);

    % Create labels to display slider values
    label1 = uilabel(fig, 'Position', [120 350 100 22], 'Text', 'Ch1 (0, 0)');
    label2 = uilabel(fig, 'Position', [380 350 100 22], 'Text', 'Ch2 (0, 0)');

    % Create the 'Done' button
    btnDone = uibutton(fig, 'Text', 'Done', 'Position', [230 20 100 22],...
        'ButtonPushedFcn', @doneButtonCallback);

    % Create buttons to toggle control
    btnCh1 = uibutton(fig, 'Text', 'Control UV', 'Position', [150 220 100 22]); 
    btnCh2 = uibutton(fig, 'Text', 'Control Green', 'Position', [400 220 100 22]);

    selectedChannel = 1; % Initialize to Channel 1 as default


    streamU = RandStream('mrg32k3a','seed',IN.NoiseUniqueSeed);
    nRows = ceil(OLED.height/(IN.NoiseGridSize/OLED.pixelSize));          %number of rows of the stimulus array
    nCols = ceil(OLED.width/(IN.NoiseGridSize/OLED.pixelSize));            %number of columns in stimulus array
    Stim_int = randi(streamU, [0, 1], nRows, nCols);
    zero_stim_int = zeros(size(Stim_int));
    % UV
    cgloadarray(1, nCols, nRows,[zero_stim_int(:),Stim_int(:),zero_stim_int(:)], OLED.width, OLED.height);
    % Green
    cgloadarray(2, nCols, nRows,[zero_stim_int(:),zero_stim_int(:),Stim_int(:)], OLED.width, OLED.height);

    cgdrawsprite(1,0,0)
    cgdrawsprite(2,0,0)
    cgtrncol(2, 'n');
    cgflip
    % Callback function for slider value changes
    function updateSliderLabels(src, ~)
        label1.Text = sprintf('UV (%d, %d)', round(hSlider1.Value), round(vSlider1.Value));
        label2.Text = sprintf('Green (%d, %d)', round(hSlider2.Value), round(vSlider2.Value));
        
         % Print the updated slider values
        x1 = round(hSlider1.Value);
        y1 = round(vSlider1.Value);
        x2 = round(hSlider2.Value);
        y2 = round(vSlider2.Value);

        % UV
        cgdrawsprite(1,x1,y1)
        cgdrawsprite(2,x2,y2)
        cgtrncol(2, 'n');
        cgflip
       
    end

    % Attach callback functions to buttons
    btnCh1.ButtonPushedFcn = @controlCh1;
    btnCh2.ButtonPushedFcn = @controlCh2;

    % Attach key press callback to the figure
    fig.KeyPressFcn = @keyPressCallback;

    % Callback function for 'Done' button
    function doneButtonCallback(~, ~)
        % Retrieve slider values before closing
        x1 = round(hSlider1.Value);
        y1 = round(vSlider1.Value);
        x2 = round(hSlider2.Value);
        y2 = round(vSlider2.Value);

        % Assign the values to varargout
        varargout{1} = x1;
        varargout{2} = y1;
        varargout{3} = x2;
        varargout{4} = y2;

        % Close the GUI
        close(fig);
    end

    % Attach callback functions to sliders
    hSlider1.ValueChangedFcn = @updateSliderLabels;
    vSlider1.ValueChangedFcn = @updateSliderLabels;
    hSlider2.ValueChangedFcn = @updateSliderLabels;
    vSlider2.ValueChangedFcn = @updateSliderLabels;

    % Button callback functions to toggle control
    function controlCh1(~, ~)
        selectedChannel = 1;
        figure(fig);
        fprintf('control UV \n')
    end

    function controlCh2(~, ~)
        selectedChannel = 2;
        figure(fig);
        fprintf('control Green \n')
    end

    % Callback function for key presses
    function keyPressCallback(~, event)
        if selectedChannel == 1  % Check if Channel 1 is selected
            switch event.Key
                case 'leftarrow'
                    hSlider1.Value = max(hSlider1.Value - 1, hSlider1.Limits(1));
                case 'rightarrow'
                    hSlider1.Value = min(hSlider1.Value + 1, hSlider1.Limits(2));
                case 'uparrow'
                    vSlider1.Value = min(vSlider1.Value + 1, vSlider1.Limits(2));
                case 'downarrow'
                    vSlider1.Value = max(vSlider1.Value - 1, vSlider1.Limits(1));
            end
        elseif selectedChannel == 2 % Check if Channel 2 is selected
            switch event.Key
                case 'leftarrow'
                    hSlider2.Value = max(hSlider2.Value - 1, hSlider2.Limits(1));
                case 'rightarrow'
                    hSlider2.Value = min(hSlider2.Value + 1, hSlider2.Limits(2));
                case 'uparrow'
                    vSlider2.Value = min(vSlider2.Value + 1, vSlider2.Limits(2));
                case 'downarrow'
                    vSlider2.Value = max(vSlider2.Value - 1, vSlider2.Limits(1));
            end
        end
        fprintf('keypressed... channel %d \n', selectedChannel)
        updateSliderLabels(); % Update labels after key press
    end

    % Bring the figure to focus
    figure(fig);

    uiwait(fig);
end
