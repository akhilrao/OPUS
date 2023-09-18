function collision_probability = calculateCollisionProbability_MOCAT4S(x0, VAR, location)
% Calculate collision probability facing a single fringe satellite in a particular location based on the state provided.
% Inputs: state vector (x0: full OUT struct from MOCAT), physical parameters (VAR: full VAR struct from MOCAT), location (column index for state struct)
% Outputs: collision probability (scalar)

k = location; % Shell index

% Define values for populations in each shell k. Always take the last element; x0 should have dimensions 1xNshells, but we want to be sure.
S = x0.S(end,k);
D = x0.D(end,k);
N = x0.N(end,k);
Su = x0.Su(end,k);

% use different probability of collision and generated number of
% fragments according to the species colliding
%% Calculate collision probabilities phi_ij
%% phi_ij is probability of collision between species i and j
%% Calculated using kinetic theory as phi_ij = pi * v_r * (r_i + r_j)^2 / V(k)
%% where v_r is relative velocity, r_i and r_j are radii, V(k) is shell volume
phi_SS = VAR.phi(1,1)/VAR.V(k);
phi_SD = VAR.phi(1,2)/VAR.V(k);
phi_SN = VAR.phi(1,3)/VAR.V(k);
phi_SSu = VAR.phi(1,4)/VAR.V(k);
phi_DD = VAR.phi(2,2)/VAR.V(k);
phi_DN = VAR.phi(2,3)/VAR.V(k);
phi_DSu = VAR.phi(2,4)/VAR.V(k);
phi_NN = VAR.phi(3,3)/VAR.V(k);
phi_NSu = VAR.phi(3,4)/VAR.V(k);
phi_SuSu = VAR.phi(4,4)/VAR.V(k);

A2_Su = -(VAR.delta+VAR.alpha)*phi_NSu; % Collisions with debris
A3_Su = -(VAR.delta+VAR.alpha)*phi_DSu; % Collisions with derelicts
A4_Su = -VAR.alpha_active*phi_SuSu; % Collisions between unslotted satellites
A5_Su = -VAR.alpha_active*phi_SSu*Su/(Su+S); % Collisions with slotted satellites

% Calculate rate of collisions for unslotted satellites in shell k
% Units: #/second
collision_rate = A2_Su*Su*N + A3_Su*Su*D + A4_Su*Su^2 + A5_Su*S*Su;

% Convert rate to probability. 
% Units: %
% Negative exponential form.
%% Notes: Since the collision rate coefficients already have negative signs in front, we don't put them down here
%% This form is used in Letizia et al (2017) and Rao and Rondina (2023). It reflects the probability of an individual satellite being destroyed in a collision, using a pigeonhold principle calculation given collision_rate many events will occur.
collision_probability = 1 - exp(collision_rate);
% Linear form.
%% Notes: Since the collision rate coefficients have negative signs in front, we need another here to make them positive
%% This form reflects the expected share of satellites that will be lost in collisions at a point in time.
% collision_probability = -collision_rate/Su;

end
