function [ydot] = DebrisDynamicsSTM_GM(t,y,GMPHD_params)
% Function to propagate debris dynamics using state transition matrix (STM) and GMPHD parameters.
% Inputs:
%   t: time vector of nodes to compute values at
%   y: state vector with Gaussians
%   GMPHD_params: struct containing the following parameters:
%       rho0: atmospheric density
%       H: scale height
%       r0: radius of Earth
%       Cd: drag coefficient
%       rhodeb: density of debris
%       mu: gravitational parameter of Earth

% Extract parameters from GMPHD_params
rho0=GMPHD_params.rho0;
H=GMPHD_params.H;
r0=GMPHD_params.r0;
Cd=GMPHD_params.Cd;
rhodeb=GMPHD_params.rhodeb;
mu=GMPHD_params.mu;

nstatesper=2+2^2;

nparam=length(y)/nstatesper;

ydot=zeros(1,length(y));

for i=1:nparam

    ii=i-1;

    yi=y((nstatesper*ii+1):(nstatesper*ii+nstatesper));

    ydot_i = DebrisDynamicsSTM(t,yi,rho0,H,r0,Cd,rhodeb,mu);

    ydot((nstatesper*ii+1):(nstatesper*ii+nstatesper)) = ydot_i;

end

ydot=ydot';

end