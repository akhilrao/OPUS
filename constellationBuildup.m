function lam_i = constellationBuildup(location_index, final_size, linear_rate, Si)
    % constellationBuildup sets the launch rate for a given constellation at a given location.
    %
    % INPUTS:
    %   location_index: Integer, the index of the location where the constellation is to be built.
    %   final_size: Integer, the final size of the constellation.
    %   linear_rate: Integer, the linear rate of buildup for the constellation.
    %   Si: Vector, the current sizes of all locations.
    %
    % OUTPUT:
    %   lam_i: Scalar, the modified launch rate for the specified location index.
    %
    % EXAMPLE USAGE:
    %   lam_i = constellationBuildup(11, 3000, 1500, Si);

    % Compute the current size of the constellation at the given location
    current_size = Si(location_index);

    % Compute the remaining size to be built
    remaining_size = max(final_size - current_size, 0);

    % Compute the modified launch rate for the given location
    lam_i = min(linear_rate, remaining_size);
end

