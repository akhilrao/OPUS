function [J] = DebrisJacobian(rdeb,adeb,rho0,H,r0,Cd,rhodeb,mu)
%Calculate Jacobian for orbital debris in circular orbits with drag, all
% units in km, kg, and s

drhoda=-(rho0/H)*exp(-(adeb-r0)/H); % density gradient at altitude a

rho=rho0*exp(-(adeb-r0)/H); % density at altitude a

J=[-4*Cd*sqrt(adeb*mu)/(3*rdeb*rhodeb)*drhoda-rho*2*Cd/(3*rdeb*rhodeb)*sqrt(mu/adeb), rho*4*Cd*sqrt(adeb*mu)/(3*rdeb^2*rhodeb);...
    0,0]; % Jacobian matrix
end