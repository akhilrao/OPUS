function [excess_returns] = excessReturnCalculator(launches, launch_mask, propagator, propagator_params, state_matrix, revenue_model, econ_params, constellation_launches)
% excessReturnCalculator - Calculate excess returns to open-access firms in each location
% Inputs:
%   launches        - Vector of open-access launch rates at each location (dimension 1 x propagator_params.N_shell)
%   launch_mask     - Vector of 1s and 0s indicating whether a location is accessible (1) or not (0) (dimension 1 x propagator_params.N_shell)
%   propagator      - Function handle for the propagator (default: MOCAT4S)
%   propagator_params - Parameters for the propagator (default: VAR)
%   state_matrix    - Object containing variables for S, Su, N, D in all shells
%   revenue_model   - Revenue model (default: 'linear')
%   econ_params  - Associated parameters for the revenue model (default: econ_params)
%   constellation_launches - Number of launches by the constellation(s) (default: lam(:,2))
% Outputs:
%   excess_returns  - Vector of excess returns for each location (dimension 1 x propagator_params.N_shell)
%
% Author: Akhil Rao
% Last updated: September 2023

%%%%%
% SET DEFAULTS

propagator_string = propagator;

% Assign function values for propagator and probability model
%% MOCAT
if strcmp(propagator_string, 'MOCAT')
    propagator = @MOCAT4S;
    probability_model = @calculateCollisionProbability_MOCAT4S;
    locations = linspace(1,40,40); % Create vector of altitude/location nodes
    tspan = linspace(0,1,2); % Create linearly-spaced vector from 0:1 with 2 points. This is the within-year time-stepper for integration.
end
%% GMPHD
if strcmp(propagator_string, 'GMPHD')
    propagator = @PropagateGM_timehist_collisions;
    probability_model = @calculateCollisionRate_GMPHD;
    t = propagator_params.t;
    GMPHD_params = propagator_params.GMPHD_params;
    shellstruct = propagator_params.shellstruct;
    GM_coll = propagator_params.GM_coll;
    i_collisions = propagator_params.i_collisions;
    GM1 = state_matrix;
end

%% Default values in case no propagator, probability model, or revenue model is provided
if isempty(propagator_string)
    propagator = @MOCAT4S; % Default propagator
    probability_model = @calculateCollisionProbability_MOCAT4S; % Default collision probability model
end
if isempty(revenue_model)
    revenue_model = 'linear'; % Default revenue model
end

%%%%%
% CALCULATE EXCESS RETURNS

locations = linspace(1,40,40); % Create vector of altitude/location nodes. Assumption is both MOCAT4S and GMPHD use 40 nodes.

% Create object to hold launch rates. Using same format as MOCAT4S setup.
lam(:,1) = constellation_launches;
lam(:,2) = launches;

% Propagate state
if strcmp(propagator_string, 'MOCAT')
    state_next_path = propagator(tspan, state_matrix, lam, propagator_params);
    state_next = state_next_path(end,:);
    % Calculate collision probability at all locations
    for i = 1:propagator_params.N_shell
        collision_probability(i) = probability_model(state_next, propagator_params, i);
    end
    % Calculate revenue using the provided revenue model and parameters. Mask is applied inside to revenues.
    rate_of_return = fringeRateOfReturn("linear", econ_params, state_next, locations, launch_mask, "MOCAT");
end

if strcmp(propagator_string, 'GMPHD')
    n_fringe_sats = zeros(40);
    dt = (t(i_collisions(2)) - t(i_collisions(1))); % this assumes only propagating for one time period
    % Assign launches guess vector to shellstruct launches slot
    [t,GM1,shellstruct] = propagator(t, GM1, GMPHD_params, shellstruct, i_collisions, GM_coll);
    GM1_last = GM1(end,:); % Get last element of the state to use for calculating collision rate
    % Calculate collision probability at all locations
    for i = 1:40
        collision_probability(i) = probability_model(GM1_last, GMPHD_params, shellstruct(i), 'fringe', 1, dt); % probability for excess return calculation

        collision_rate(i) = probability_model(GM1_last, GMPHD_params, shellstruct(i), 'fringe', 0, dt); % rate for law of motion
        ncoll_fringe(i)=min([collision_rate(i), shellstruct(i).nfringe(end)]*(1 - 1/GMPHD_params.satLifetime)); % Keep it within physical limits

        % shellstruct(i).nfringe(end+1) = shellstruct(i).nfringe(end)*(1 - 1/GMPHD_params.satLifetime) - ncoll_fringe(i) + launches(i); % Update number of fringe satellites

        % n_fringe_sats(i) = shellstruct(i).nfringe(end);

        n_fringe_sats(i) = shellstruct(i).nfringe(end)*(1 - 1/GMPHD_params.satLifetime) - ncoll_fringe(i) + launches(i); 
    end
    % Calculate revenue using the provided revenue model and parameters. Mask is applied inside to revenues.
    rate_of_return = fringeRateOfReturn("linear", econ_params, n_fringe_sats, locations, launch_mask, "GMPHD");
end

% Calculate the excess returns
excess_returns = 100*(rate_of_return - collision_probability*(1 + econ_params.tax));