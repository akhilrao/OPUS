function [GM_init]=prepareInitialGMs(GMPHD_params, nyears, ntime,dicoll)
% Function to prepare time step nodes and initial state object with Gaussians, GM, for GMPHD propagator. Takes in GMPHD parameters and time horizon parameters.
% Inputs:
%   GMPHD_params: GMPHD parameters
%   nyears: number of years to propagate
%   ntime: number of time steps to propagate
% Outputs:
%   GM_init: a struct containing the following objects:
%       t: time vector of nodes to compute values at
%       GM: initial state object
%       i_collisions: time steps to compute collisions at
%       GM_coll: number of debris per collision

    i_collisions=[1,dicoll:dicoll:ntime]; % Time steps for collision calculations

    % Set parameters from params element
    nr=GMPHD_params.nr; % number of radius nodes
    na=GMPHD_params.na; % number of altitude nodes
    nt=GMPHD_params.nt; % product of altitude and radius nodes
    rE=GMPHD_params.rE; % radius of Earth

    Cd=GMPHD_params.Cd; % unitless
    rhodeb=GMPHD_params.rhodeb; % kg/m^3

    mu=GMPHD_params.mu; % m^3/s^2
    H=GMPHD_params.H; % scale height
    rho0=GMPHD_params.rho0; % kg/m^3

    aarr=GMPHD_params.aarr; % array of debris altitudes
    rarr=GMPHD_params.rarr; % array of debris radii
    P0=GMPHD_params.P0; % covariance matrix for debris
    r0=GMPHD_params.r0;

    % Initialize state object
    GM0.state=[];
    GM0.P=[];
    GM0.c=[];

    % Debris per collision
    GM_coll(1)=GMPHD_params.GM_coll(1);
    GM_coll(2)=GMPHD_params.GM_coll(2);
    GM_coll(3)=GMPHD_params.GM_coll(3);
    GM_coll(4)=GMPHD_params.GM_coll(4);

    % Initialize GM for debris
    GM(nt)=GM0;
    ii=1;
    for i=1:na
        for j=1:nr
            GM(ii).state=[aarr(i),rarr(j)];
            GM(ii).P=P0;
            GM(ii).c=5000;
        ii=ii+1;
        end
    end

    t=linspace(0,nyears*365.25*24*60*60,ntime);

    GM_init.t = t;
    GM_init.GM = GM;
    GM_init.i_collisions = i_collisions;
    GM_init.GM_coll = GM_coll;

end