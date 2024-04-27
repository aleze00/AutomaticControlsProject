%% Fix 'OutDataTypeStr' for 'myVars' bus signals
% Define the path to your Simulink model
model_path = 'NewMODEL';
DataType = 'Bus: slBus1'
% Open the Simulink model
open_system(model_path);

% Get all blocks in the model
blocks = find_system(model_path, 'SearchDepth', 1, 'LookUnderMasks', 'all', 'Type', 'Block');

% Iterate through each block
for i = 1:length(blocks)
    block = blocks{i};
    % Check if the block is a Constant block
    if strcmp(get_param(block, 'BlockType'), 'Constant')
        % Check if the variable name is 'myAs'
        variableName = get_param(block, 'Value');
        if strcmp(variableName, 'myAs')
            % Set the 'OutDataTypeStr' parameter to 'Bus: slBus1'
            set_param(block, 'OutDataTypeStr', DataType);
        end
    end
end

% Save the modified model
save_system(model_path);

% Close the Simulink model
close_system(model_path);


