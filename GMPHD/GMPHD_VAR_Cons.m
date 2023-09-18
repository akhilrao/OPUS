function [GMPHD_params] = GMPHD_VAR_Cons()

% Constants
GMPHD_params.rE=6378; % Earth radius [km]
GMPHD_params.Cd=2.2; % Drag coefficient
GMPHD_params.rhodeb=1e11;
GMPHD_params.H=88.667; % Atmospheric density scale height [km]
GMPHD_params.rho0=3.614e-5; % Atmospheric density at h=0 [kg/km^3]
GMPHD_params.mu=398600; % Earth gravitational parameter [km^3/s^2]
GMPHD_params.rdebmin=0;
GMPHD_params.rdebmax=50;
GMPHD_params.satLifetime=5; % Satellite lifetime [years]

% Parameters needed for GMPHD discretization
%% Update the nr and na parameters for computational fidelity/tractability
GMPHD_params.nr=8; % number of debris object radius nodes
GMPHD_params.na=8; % number of altitude nodes
GMPHD_params.nt=GMPHD_params.nr*GMPHD_params.na; % product for grid
GMPHD_params.aarr=linspace(GMPHD_params.rE+300,GMPHD_params.rE+950,GMPHD_params.na); % array of altitudes / semi-major axes
GMPHD_params.rarr=linspace(0.01/1000,1/1000,GMPHD_params.nr); % array of debris radii

% Parameters needed for GMPHD distributions
GMPHD_params.P0=[50, 0; 0, 0.5/1000].^2; % Debris covariance

% Parameters needed for GMPHD propagation
%% Debris per collision -- breakup model. All collisions assumed to be catastrophic.
GMPHD_params.GM_coll(1).state=[0,0.1/1000]; % first element is altitude (SMA), second is radius of debris -- this is very small stuff
GMPHD_params.GM_coll(1).P=[20^2,0;0,(0.04/1000)^2]; % covariance of debris distribution, measures how far the debris spreads following a collision. 20 is set for 100 km shells, go smaller to keep them inside 35 km shells.
GMPHD_params.GM_coll(1).c=40; % Cardinality of debris distribution, what the random finite set Gaussian thing integrates to

GMPHD_params.GM_coll(2).state=[0,0.25/1000];
GMPHD_params.GM_coll(2).P=[20^2,0;0,(0.08/1000)^2];
GMPHD_params.GM_coll(2).c=30;

GMPHD_params.GM_coll(3).state=[0,0.4/1000];
GMPHD_params.GM_coll(3).P=[20^2,0;0,(0.1/1000)^2];
GMPHD_params.GM_coll(3).c=20;

GMPHD_params.GM_coll(4).state=[0,0.7/1000];
GMPHD_params.GM_coll(4).P=[20^2,0;0,(0.1/1000)^2];
GMPHD_params.GM_coll(4).c=10;
