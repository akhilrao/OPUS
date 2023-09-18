function [VAR, econ_params, constellation_params] = modifyParameters(VAR, econ_params, constellation_params, parameter_file)
    % modifyParameters - Modifies the parameters VAR and econ_params based on a CSV file
    %
    % Inputs:
    %   VAR - Struct containing physical parameters
    %   econ_params - Struct containing economic parameters
    %   constellation_params - Struct containing constellation parameters
    %   parameter_file - CSV file containing parameter modifications
    %
    % Outputs:
    %   VAR - Modified physical parameters
    %   econ_params - Modified economic parameters
    %   constellation_params - Modified constellation parameters
    %
    % NOTES: Code currently assumes that the parameter_value column contains numerical values. If any of the parameters are meant to be strings or other non-numeric types, additional handling will be required.

    % Check if parameter_file is not "benchmark"
    if ~strcmp(parameter_file, 'benchmark')
        % Read the CSV file
        parameters = readtable(parameter_file);

        % Loop through each row of the parameters table
        for i = 1:height(parameters)
            parameter_type = parameters.parameter_type{i};
            parameter_name = parameters.parameter_name{i};
            parameter_value = parameters.parameter_value(i);

            % Modify the value based on parameter_type
            switch parameter_type
                case 'MOCAT'
                    % If the field exists in VAR, update its value
                    if isfield(VAR, parameter_name)
                        VAR.(parameter_name) = parameter_value;
                    end
                case 'econ'
                    % If the field exists in econ_params, update its value
                    if isfield(econ_params, parameter_name)
                        econ_params.(parameter_name) = parameter_value;
                    end
                otherwise
                    warning('Unknown parameter_type: %s', parameter_type);
            end
        end
    end
end
