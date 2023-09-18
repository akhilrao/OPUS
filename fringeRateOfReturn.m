function rateOfReturn = fringeRateOfReturn(revenue_model, econ_params, state_next, target_location, mask, model_type)
% Fringe satellite rate of return function. 
% Inputs: model type (revenue_model), params (econ_params), state (state_next), target shell (target_location).
%   revenue_model: string ,"linear" or "CES"
%   econ_params: struct, contains parameters for revenue model
%   state_next: struct, contains state variables for next time step
%   target_location: scalar, shell index of target location. Can accept a scalar or a vector. If a vector, will return a vector of rate of return values for each index.
%   mask: binary vector, mask of which shells are accessible
%   model_type: string, "MOCAT" or "GMPHD"
% Outputs: rate of return on cost of deploying to target_location. May output a vector of dimension [target_location] if target_location is a vector rather than a scalar
% Output units: [%/satellite/year]
% Notes: Start with a linear form with a common coefficient for all shells. Divide by shell-specific cost coefficient to put this into "rate" units. All choice logic for launchers deploying to specific shells is embedded in here.
%
% Author: Akhil Rao
% Last updated: September 2023

% Calculate net rate of return based on the given revenue model and parameters
if strcmp(revenue_model, 'linear')
    % Linear revenue model: intercept - coef * total_number_of_fringe_satellites
    intercept = econ_params.intercept;
    coef = econ_params.coef;

    if strcmp(model_type, 'MOCAT')
        fringe_total = state_next.Su(end,:); % Total number of fringe satellites across shells
    end
    if strcmp(model_type, 'GMPHD')
        fringe_total = state_next; % Total number of fringe satellites across shells
    end

    revenue = intercept - coef * sum(fringe_total); % Revenue is a linear decreasing function of the total number of fringe satellites, common coefficients across shells
    revenue = revenue .* mask(target_location); % Only return revenue for shells that are accessible 

    cost = econ_params.cost(target_location); % Pull cost schedule from the params struct. Should be a vector of length N_shells

    discount_rate = econ_params.discountRate; % Discount rate, scalar

    depreciation_rate = 1/econ_params.satLifetime; % 1 / lifetime, scalar

    rateOfReturn = revenue ./ cost - discount_rate - depreciation_rate; % Equilibrium expression for rate of return.
else
    % Other revenue models can be implemented here
    rateOfReturn = 0; % Placeholder value
end