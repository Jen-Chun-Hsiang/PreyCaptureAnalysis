function run(obj)
    % Main method to execute the protocol
    lightCrafter = LightCrafterDevice(yourSettings); % Configure with your settings
    movingBar = MovingBar(); % Initialize your moving bar stimulus
    
    obj.prepareRun();
    for i = 1:obj.totalNumEpochs
        disp(['Starting epoch ' num2str(i) ' of ' num2str(obj.totalNumEpochs)]);
        obj.prepareEpoch();
        % Here, integrate the moving bar stimulus presentation
        movingBar.run(); % This will display the moving bar at the configured settings
        pause(obj.stimTime / 1000); % Adjust timing as necessary
        disp(['Completed epoch ' num2str(i)]);
    end
    obj.completeRun();
end
