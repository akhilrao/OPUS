function [t,GMP,shellstruct] = PropagateGM_timehist_collisions(t, GM0, GMPHD_params, shellstruct, i_collisions, GMr_coll)
% Function to propagate state of the debris environment using Gaussian Mixture Particle Hypothesis Density (GMPHD) filter.
% Inputs: 
%   t: time vector of nodes to compute values at
%   GM0: Gaussian mixtures of debris state at t(1)
%   GMPHD_params: struct containing the following parameters:
%       rho0: atmospheric density
%       H: scale height
%       r0: radius of Earth
%       Cd: drag coefficient
%       rhodeb: density of debris
%       mu: gravitational parameter of Earth
%       rdebmin: minimum radius of debris
%       rdebmax: maximum radius of debris
%   shellstruct: structure containing shell parameters. It has the following fields: rmin, rmax, nfringe, nconst, rfringe, rconst, ffailconst, ffailfringe, V, launch rates for fringe and const. See buildShellStruct.m for more details.
%   i_collisions: vector of indices of collision timesteps. Begins with 1 even though no collisions occur at the first timestep.
%   GMr_coll: Gaussian mixture of debris state at collision
% Outputs:
%   t: time vector of nodes that values were computed at (same as input)
%   GMP: Gaussian mixtures of debris state at t
%   shellstruct: structure containing shell parameters

I=eye(2); % Identity matrix

GMP(1,length(GM0))=GM0(1); % Initialize GMP

% Loop over collision timestep intervals
for kk=2:length(i_collisions)

    y=[]; % Initialize state vector

    % Populate state vector
    for ii=1:length(GM0)
        y=[y, GM0(ii).state, reshape(I,1,4)];
    end

    % Propagate debris dynamics
    % tic;
    % fprintf('Now starting GMPHD propagation...\n');
    warning('off', 'all');
    tolerance=1e-14;
    options = odeset('RelTol',tolerance,'AbsTol',tolerance);
    [~,yfi] = ode45(@DebrisDynamicsSTM_GM,t(i_collisions(kk-1):i_collisions(kk)), y, options, GMPHD_params);
    warning('on', 'all');
    % fprintf('GMPHD propagation complete. Elapsed time: %f seconds.\n', toc);

    % Update Gaussian mixture within each collision timestep
    for j=i_collisions(kk-1):i_collisions(kk)

        yf=squeeze(yfi(j-i_collisions(kk-1)+1,:)); % Extract state vector

        % Loop over each Gaussian to update state and covariance
        for i=1:length(GM0)
            ii=i-1; % Index for previous timestep
            GMP(j,i).state=yf(ii*6+1:ii*6+2); % Update state
            phi=reshape(yf(ii*6+3:ii*6+6),2,2); % Extract state transition matrix
            P=GM0(i).P; % Extract covariance matrix
            P=phi*P*phi'; % Propagate covariance matrix
            GMP(j,i).P=P; % Update covariance matrix
            GMP(j,i).c=GM0(i).c; % Copy weight

        end
    end

    GM1(1,1)=GM0(1); % Initialize GM1

    c=1; % Counter

    % Remove Gaussians with infinite covariance
    for i=1:length(GM0)
        % Check if covariance is finite
        if isfinite(GMP(end,i).state(1))
                GM1(c).state=GMP(end,i).state; % Copy state
                GM1(c).c=GMP(end,i).c; % Copy weight
                GM1(c).P=GMP(end,i).P; % Copy covariance
            c=c+1; % Increment counter
        end
    end

    GM0=GM1; % Update GM0
    clear GM1; % Clear GM1

    dt=(t(i_collisions(kk))-t(i_collisions(kk-1))); % Time step

    % Calculate collisions in each shell
    for i=1:length(shellstruct)

        cdot_total_fringe = calculateCollisionRate_GMPHD(GM0, GMPHD_params, shellstruct(i), 'fringe', 0, dt);
        ncoll_fringe=min([cdot_total_fringe, shellstruct(i).nfringe(end)]*(1 - 1/GMPHD_params.satLifetime)); % Keep it within physical limits

        cdot_total_const = calculateCollisionRate_GMPHD(GM0, GMPHD_params, shellstruct(i), 'const', 0, dt);
        ncoll_const=min([cdot_total_const, shellstruct(i).nconst(end)]*(1 - 1/GMPHD_params.satLifetime)); % Keep it within physical limits

        shellstruct(i).ncoll(end+1) = ncoll_fringe+ncoll_const; % Store total number of collisions
        shellstruct(i).nconst(end+1) = shellstruct(i).nconst(end)*(1 - 1/GMPHD_params.satLifetime) - ncoll_const + shellstruct(i).launch_rate_const; % Update number of constellation satellites
        shellstruct(i).nfringe(end+1) = shellstruct(i).nfringe(end)*(1 - 1/GMPHD_params.satLifetime) - ncoll_fringe + shellstruct(i).launch_rate_fringe; % Update number of fringe satellites

        GMColl=GMr_coll; % Copy Gaussian mixture

        % Update Gaussian mixture
        for j=1:length(GMColl)
            GMColl(j).state(1) = (shellstruct(i).rmin + shellstruct(i).rmax) / 2; % Set radius to middle of shell
            GMColl(j).c = GMColl(j).c * (ncoll_fringe + ncoll_const); % Update weight
        end

        GM0((end+1):(end+length(GMColl)))=GMColl; % Add to GM0

    end


end


end