function scenarioNamer(parameter_file, model_type)
    % scenarioNamer - Generates a unique and memorable name for a scenario based on the given parameters
    % Inputs:
    %   parameter_file - CSV file containing parameter modifications. If "benchmark," default parameters are used.
    % Outputs:
    %   unique_name - Unique and memorable name for the scenario

    %% Add path to MOCAT4S and GMPHD filter files
    addpath('MOCAT4S')
    addpath('GMPHD')

    % Set up physical model parameters from MOCAT
    VAR = MOCAT4S_VAR_Cons();

    % Set up GMPHD parameters
    GMPHD_params = GMPHD_VAR_Cons();

    % Set up economic model parameters
    econ_params = set_econ_parameters(VAR);

    %%%% CONSTELLATION BLOCK
    % Set up constellation parameters. Read in from constellation-parameters.csv in parsets folder. Change this to run different constellation scenarios.
    constellation_params = constellation_parameters('scenarios/parsets/constellation-parameters.csv');
    %%%% END CONSTELLATION BLOCK

    %%%% LOCATION MASKING BLOCK
    % TODO: Make this a switch based on string input. Definition also appears in iam_solver.m
    % Create mask(s): 1 if location is under open access, 0 if not. Naming convention is [where it's 1]-mask
    % unmask: 1 for all shells
    launch_mask = ones(1,40);
    %%%% END LOCATION MASKING BLOCK

    % Modify parameters if a specific parameter file is provided (not "benchmark")
    if ~strcmp(parameter_file, 'benchmark')
        [VAR, econ_params, constellation_params] = modifyParameters(VAR, econ_params, constellation_params, parameter_file);
    end

    % Serialize the VAR and econ_params structs into JSON strings
    jsonStr_VAR = jsonencode(VAR);
    jsonStr_GMPHD_params = jsonencode(GMPHD_params);
    jsonStr_econ_params = jsonencode(econ_params);
    jsonStr_constellation_params = jsonencode(constellation_params);
    jsonStr_launch_mask = jsonencode(launch_mask);

    % Concatenate the JSON strings
    combinedJsonStr = [jsonStr_VAR, jsonStr_GMPHD_params, jsonStr_econ_params, jsonStr_constellation_params, jsonStr_launch_mask];

    % Generate a hash value from the combined JSON string
    md = java.security.MessageDigest.getInstance('SHA-256');
    hash = md.digest(uint8(combinedJsonStr));
    hash_value = sprintf('%02x', hash);

    % Generate a unique name
    unique_name = generateUniqueName(hash_value);
    unique_name

    % Write the unique name to a text CSV file with the unique name as filename and scenario table as contents. The file is saved in the "scenarios" folder under a sub-folder with the same unique name. The sub-folder must be created first.
    mkdir(['scenarios/', unique_name]);
    % writetable(scenario_table, ['scenarios/', unique_name, '/', unique_name, '--parameters.csv']);
    % Write out a text file named "scenario_name.txt" with the unique name
    fileID = fopen('scenario_name.txt', 'w');
    fprintf(fileID, '%s', unique_name);
    fclose(fileID);

end
