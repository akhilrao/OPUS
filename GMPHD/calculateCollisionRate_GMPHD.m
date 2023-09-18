function [cdot_total]=calculateCollisionRate_GMPHD(GM0, GMPHD_params, shellstruct, sat_type, probability, dt)
    % Function to wrap calcCdGaussianStruct and calculate number of collisions and collision rate for fringe or constellation satellites
    % Inputs: 
    %   GM0: Gaussian mixtures of debris state
    %   GMPHD_params: struct containing the following parameters:
    %       rmin: minimum radius of shell
    %       rmax: maximum radius of shell
    %       rdebmin: minimum radius of debris
    %       rdebmax: maximum radius of debris
    %       V: volume of shell
    %       r: radius of shell
    %       n: number of satellites in shell
    %       ffail: failure probability
    %   shellstruct: structure containing shell parameters
    %   sat_type: 'fringe' or 'const'
    %   probability: flag for whether to use negative exponential model to convert rate to probability or not
    %   dt: length of time interval
    % Outputs:
    %   cdot_total: total number of collisions
    %   cdot_ind: number of collisions for each Gaussian

    % Extract parameters from GMPHD_params
    rho0=GMPHD_params.rho0;
    H=GMPHD_params.H;
    r0=GMPHD_params.r0;
    Cd=GMPHD_params.Cd;
    rhodeb=GMPHD_params.rhodeb;
    mu=GMPHD_params.mu;
    rdebmin=GMPHD_params.rdebmin;
    rdebmax=GMPHD_params.rdebmax;

    % Number of collisions in fringe
    if strcmp(sat_type,'fringe')
        [cdot_total,~]=calcCdGaussianStruct(GM0, shellstruct.rmin, shellstruct.rmax, rdebmin,rdebmax, shellstruct.V, shellstruct.rfringe, shellstruct.nfringe(end), shellstruct.ffailfringe); % Calculate from Gaussians
    end

    % Number of collisions in constellation
    if strcmp(sat_type,'const')
        [cdot_total,~]=calcCdGaussianStruct(GM0,shellstruct.rmin, shellstruct.rmax, rdebmin,rdebmax, shellstruct.V,shellstruct.rconst,shellstruct.nconst(end), shellstruct.ffailconst); % Calculate from Gaussians
    end
    
    cdot_total = cdot_total*dt; %Note: don't want to apply the physical limits here, since the rate is a continuous-time object, while the negative exponential used in the probability calculation is a discrete-time object that only goes to 1 at infinity.

    % Convert rate to probability
    if probability
        cdot_total=1-exp(-cdot_total); % Negative exponential model of probability of collision within time interval dt
    end

end