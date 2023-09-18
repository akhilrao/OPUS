function [launch_rate] = openAccessSolver(launch_rate_input, launch_mask, propagator, propagator_params, state_matrix, revenue_model, econ_params, constellation_launches, n_workers)
% openAccessSolver - Calculate open-access launch rate
% Inputs:
%   launch_rate_input - Initial guess or solver mask for open-access launch rates (dimension 1 x propagator_params.N_shell)
%   launch_mask     - Vector of 1s and 0s indicating whether a location is accessible (1) or not (0) (dimension 1 x propagator_params.N_shell)
%   propagator      - Function handle for the propagator (default: MOCAT4S)
%   propagator_params - Parameters for the propagator (default: VAR)
%   state_matrix    - Object containing variables for S, Su, N, D in all shells (MOCAT) or debris Gaussians (GMPHD)
%   revenue_model   - Revenue model (default: 'linear')
%   econ_params  - Associated parameters for the revenue model (default: econ_params)
%   constellation_launches - Number of launches by the constellation(s) (default: lam(:,2))
%   n_workers       - Number of workers for parallel computing (default: 1)
% Outputs:
%   launch_rate     - Open-access launch rates (dimension 1 x propagator_params.N_shell)
%
% Author: Akhil Rao
% Last updated: September 2023

% Nested function to pass additional parameters to excessReturnCalculator
function y = excessReturn(launches)
    y = excessReturnCalculator(launches, launch_mask, propagator, propagator_params, state_matrix, revenue_model, econ_params, constellation_launches);
end

    % Initialize parallel pool iff n_workers > 1
    if n_workers > 1
        parpool('local', n_workers)
    end

    % Reassignment to allow for masking initial guess
    launch_rate_init = launch_rate_input.*launch_mask;

    lower_bound = zeros(size(launch_rate_init)); % Lower bound for launch rates

    % Solve system of equations for open-access equilibrium
    options = optimoptions('lsqnonlin', 'UseParallel', true, 'Display', 'iter', 'Algorithm', 'levenberg-marquardt'); % Display iteration information
    [launch_rate, res] = lsqnonlin(@excessReturn, launch_rate_init, lower_bound, [], options); % Solve system of equations

    p = gcp('nocreate'); % Get the current parallel pool (if it exists)
    if ~isempty(p)
        delete(p); % Shut down the parallel pool
    end

end