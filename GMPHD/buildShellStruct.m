function [shellstruct]=buildShellStruct(const_levels, fringe_levels, const_launches, fringe_launches, GMPHD_params)
    % Function to build shellstruct, called buildShellStruct. Takes in vectors of constellation and fringe satellite levels, vectors of fringe and constellation satellite launch rates. Returns shellstruct, a structure array with the following fields:
    % rmin: lower radius of shell [km]
    % rmax: upper radius of shell [km]
    % nconst: number of constellation satellites in the shell initially
    % nfringe: number of fringe satellites in the shell initially
    % ffailconst: fraction of constellation satellites that fail to avoid collisions per time step. Hard-coded to 0.01, i.e. 99% successful collision avoidance
    % ffailfringe: fraction of fringe satellites that fail to avoid collisions per time step. Hard-coded to 0.05, i.e. 95% successful collision avoidance
    % rconst: radius of constellation satellites [km]
    % rfringe: radius of fringe satellites [km]
    % launch_rate_const: constellation satellites launched per time step
    % launch_rate_fringe: fringe satellites launched per time step
    % V: relative velocity [m/s]
    % ncoll: initial number of collisions per satellite per time step
    
    for i=1:40
        shellstruct(i).rmin=GMPHD_params.rE+300+(i-1)*35; % lower radius of shell [km]
        shellstruct(i).rmax=GMPHD_params.rE+300+i*35; % upper radius of shell [km]
        % shellstruct(i).rmin=rE+200+(i-1)*35; 
        % shellstruct(i).rmax=rE+200+i*35; 
    
        shellstruct(i).nconst=const_levels(i); % number of constellation satellites in the shell initially
        shellstruct(i).nfringe=fringe_levels(i); % number of fringe satellites in the shell initially
        shellstruct(i).ffailconst=0.01; % fraction of constellation satellites that fail to avoid collisions per time step
        shellstruct(i).ffailfringe=0.05; % fraction of fringe satellites that fail to avoid collisions per time step
    
        shellstruct(i).rconst=5/1000; % radius of constellation satellites [km]
        shellstruct(i).rfringe=5/1000; % radius of fringe satellites [km]
    
        shellstruct(i).launch_rate_const=const_launches(i); % constellation satellites launched per time step
        shellstruct(i).launch_rate_fringe=fringe_launches(i); % fringe satellites launched per time step
    
        shellstruct(i).V=10; % relative velocity [m/s]
        shellstruct(i).ncoll=[0]; % initial number of collisions per satellite per time step
    end