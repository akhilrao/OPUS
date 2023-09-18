function [costFnParams] = buildCostFunction(VAR, econParams)
%%%%%
% Function to build satellite cost function using physical and economic parameters. Uses the disposal time regulation in the economic parameters struct to calculate the highest compliant altitude, then uses the delta-v required for Hohmann transfer to construct a cost function
%%%%%

%%%% BEGIN DISPOSAL COMPLIANCE ALTITUDE CALCULATION BLOCK
%% Calculates the cumulative residence time at each altitude, maps it to shell indices, then creates new variables to store the maximum altitudes compliant with the disposal time regulation from set_econ_parameters. The calculated index of the highest naturally-compliant altitude, k_star, is later appended to VAR and econParams in iam_solver.
% Initialize vectors
shell_marginal_decay_rates = zeros(VAR.N_shell, 1);
shell_marginal_residence_times = zeros(VAR.N_shell, 1);
shell_cumulative_residence_times = zeros(VAR.N_shell, 1);
% Loop over shells and compute marginal and cumulative residence times
for k=1:VAR.N_shell
    rhok = densityexp(VAR.R02(k)); % no solar flux
    rvel_current_D = -rhok*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25); % relative velocity
    shell_marginal_decay_rates(k) = -rvel_current_D/VAR.Dhl; %negative is put in to undo the negative in rvel_current_D
    shell_marginal_residence_times(k) = 1 / shell_marginal_decay_rates(k); %residence time for derelicts computed as inverse decay rate
end
shell_cumulative_residence_times = cumsum(shell_marginal_residence_times);
% Find the index of shell_cumulative_residence_times, k_star, which is the largest index such that shell_cumulative_residence_times(k_star) <= econ_params.disposalTime
indices = find(shell_cumulative_residence_times <= econParams.disposalTime);
costFnParams.k_star = max(indices); % this is the highest altitude that is compliant with the specified disposal time regulation
%%%% END DISPOSAL COMPLIANCE ALTITUDE CALCULATION BLOCK

%%% Physics-based cost function block: uses delta-v and deorbit equirements
%% Cost function is constructed as the sum of three terms: 
% 1. "lift price": price of getting 1 kg to LEO; 
% 2. "lifetime loss cost": opportunity cost of lost revenues due to deorbit cutting life short
% 3. "stationkeeping cost": financial cost of supplying the stationkeeping delta-v, monetized at econParams.delta_v_cost; 
% CHECK: work out whether this prices all the delta-v or leaves some of the safety margins unpriced. Don't want to double count the deorbit maneuver either.
%%%% 1. BEGIN DRAG DELTA-V CALCULATION
% Initialize vector
v_drag = zeros(VAR.N_shell, 1); % delta-v needed to offset drag for stationkeeping maneuvers
% Time interval over which the drag acts
delta_t = 24*3600*365.25; % Time interval in seconds
% Loop over shells
for k=1:VAR.N_shell
    rhok = densityexp(VAR.R02(k)); % drag force density in kg/m^3
    orbital_velocity_m_s = sqrt(VAR.mu / VAR.R0(k) ); % orbital velocity in m/s. Scale by 1e-3 to get km/s
    % Calculate drag force F_drag in N
    F_drag = VAR.Cd * 0.5 * rhok * orbital_velocity_m_s^2 * VAR.A(2); % Needs to have factor of 1e-6 here if using rho from rho0*exp(-(ri-r0atm)/H)
    % Calculate annual drag acceleration in km/s^2
    v_drag(k) = F_drag / VAR.mass(2) * delta_t * 1e-3;
end

%% Plot v_drag over k on log scale
% semilogy(VAR.R02(1:40), v_drag, 'LineWidth', 2)
% xlabel('Altitude [km]')
% ylabel('Drag delta-v [km/s]')
% title('Drag delta-v vs. altitude')
% grid on
% saveas(gcf, 'v_drag.png')

%%%% 1. END DRAG DELTA-V CALCULATION
%%%% 2. BEGIN LIFETIME LOSS DUE TO DEORBIT CALCULATION
% Calculate delta-v required for deorbiting
original_orbit_delta_v = zeros(VAR.N_shell, 1);
target_orbit_delta_v = zeros(VAR.N_shell, 1);
r2 = VAR.R0(costFnParams.k_star); % target altitude
for k=1:VAR.N_shell
    original_orbit_delta_v(k) = sqrt(VAR.mu / VAR.R0(k)) * (1 - sqrt(2 * r2 / (VAR.R0(k) + r2)));
    target_orbit_delta_v(k) = sqrt(VAR.mu / r2) * (sqrt(2 * VAR.R0(k) / (VAR.R0(k) + r2))-1);
end
original_orbit_delta_v = max(0, original_orbit_delta_v);
target_orbit_delta_v = max(0, target_orbit_delta_v);
total_deorbit_delta_v = original_orbit_delta_v + target_orbit_delta_v;
%% Calculate delta-v budget for mission
delta_v_budget = 1.5 * econParams.satLifetime * v_drag + 100; % safety margin factor of 50% to deal with solar activity fluctuations; another margin of 100 m/s to allow for additional/discretionary maneuvers. [delta-v] units
%% Create indicator vector for altitudes that are naturally compliant with the given deorbit regulation
naturally_compliant_vector = zeros(VAR.N_shell, 1);
naturally_compliant_vector(1:costFnParams.k_star) = 1;
%% Calculate delta-v leftover after deorbiting [delta-v] units
delta_v_after_deorbit = max(0,delta_v_budget - total_deorbit_delta_v.*(1-naturally_compliant_vector));
%% Calculate remaining lifetime after deorbit. For orbits below or at k_star, this should just be econParams.satLifetime. For orbits above k_star, this should be (delta_v_after_deorbit/delta_v_budget)*econParams.satLifetime
lifetime_after_deorbit = max(0,(delta_v_after_deorbit./delta_v_budget).*econParams.satLifetime);
%% Calculate lifetime loss due to deorbit. [%] units
lifetime_loss = (econParams.satLifetime - lifetime_after_deorbit)/econParams.satLifetime;
%%%% 3. END LIFETIME LOSS DUE TO DEORBIT CALCULATION
%%%% BEGIN COST FUNCTION COMPILATION
% 1. "lift price": default $5000/kg lift cost (based on Corrado, Cropper, Rao (2023)); 
total_lift_price = econParams.lift_price*VAR.mass(2);
% 2. "lifetime loss cost"
lifetime_loss_cost = lifetime_loss*econParams.intercept; % This is an upper bound since it uses the max possible annual revenues to monetize the loss
deorbit_maneuver_cost = total_deorbit_delta_v*econParams.delta_v_cost; % "deorbit maneuver cost" for reporting purposes
% 3. "stationkeeping cost" ("delta-v cost")
stationkeeping_cost = delta_v_budget*econParams.delta_v_cost;
% FIN: the cost of deploying a satellite to a given altitude, accounting for deorbit compliance (VAR.P is non-compliance rate)
costFnParams.cost = (total_lift_price + stationkeeping_cost + (lifetime_loss_cost)*(1 - VAR.P))';

costFnParams.total_lift_price = total_lift_price;
costFnParams.deorbit_maneuver_cost = deorbit_maneuver_cost;
costFnParams.stationkeeping_cost = stationkeeping_cost;
costFnParams.lifetime_loss_cost = lifetime_loss_cost;
costFnParams.v_drag = v_drag;
