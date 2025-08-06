function varargout  = sliderGUI()
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

    % Callback function for slider value changes
    function updateSliderLabels(src, ~)
        label1.Text = sprintf('UV (%d, %d)', round(hSlider1.Value), round(vSlider1.Value));
        label2.Text = sprintf('Green (%d, %d)', round(hSlider2.Value), round(vSlider2.Value));

         % Print the updated slider values
        x1 = round(hSlider1.Value);
        y1 = round(vSlider1.Value);
        x2 = round(hSlider2.Value);
        y2 = round(vSlider2.Value);
        fprintf('UV: (%d, %d)\n', x1, y1);
        fprintf('Green: (%d, %d)\n', x2, y2);
    end

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

    uiwait(fig);
end
