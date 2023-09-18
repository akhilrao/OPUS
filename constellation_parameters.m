function constellation_params = constellation_parameters(filename)
    % constellation_parameters - creates a struct of constellation parameters.
    %
    % INPUTS:
    %   n_constellations: Integer, number of constellations.
    %   location_indices: Vector of integers, location indices for each constellation.
    %   sizes: Vector of integers, final sizes for each constellation.
    %   rates: Vector of integers, linear rates for each constellation.
    %
    % OUTPUT:
    %   constellation_params: Struct, contains the parameters for the constellations.
    %
    % EXAMPLE USAGE:
    %   params = constellation_parameters('scenarios/parsets/constellation-parameters.csv');
    %   disp(params.location_index);  % Outputs: [11, 29]
    %   disp(params.final_size);      % Outputs: [5000, 660]
    %   disp(params.linear_rate);     % Outputs: [1500, 300]

    % filename = 'scenarios/parsets/constellation-parameters.csv'

    % Read the constellation parameters from the input file
    initial_struct = readtable(filename);
    initial_struct = table2struct(initial_struct); % Convert to struct
    %%% Constellation characteristics
    % Define location, final size, and linear rate for constellation buildup
    n_constellations = initial_struct(1).n_constellations; % Number of constellations

    % Define the parameters for each constellation
    constellation_location_indices = [];
    sizes = [];
    rates = [];
    for(i = 1:n_constellations)
        constellation_location_indices = [constellation_location_indices, initial_struct(i).location_indices];
        sizes = [sizes, initial_struct(i).target_sizes];
        rates = [rates, initial_struct(i).max_launch_rates];
    end
    
    % Initialize vectors for the constellation parameters
    location_index = zeros(1, n_constellations);
    final_size = zeros(1, n_constellations);
    linear_rate = zeros(1, n_constellations);

    % Assign values to the vectors based on the input arguments
    location_index(:) = constellation_location_indices;
    final_size(:) = sizes;
    linear_rate(:) = rates;

    % Create a struct to store the constellation parameters
    constellation_params = struct('n_constellations', n_constellations, ...
        'location_index', location_index, ...
        'final_size', final_size, ...
        'linear_rate', linear_rate);
end
