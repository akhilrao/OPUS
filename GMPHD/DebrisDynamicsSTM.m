function [ydot] = DebrisDynamicsSTM(t,y,rho0,H,r0,Cd,rhodeb,mu)

rdeb=y(2);
adeb=y(1);

rho=rho0*exp(-(adeb-r0)/H);

ydot_state=[-rho*4*Cd*sqrt(adeb*mu)/(3*rdeb*rhodeb),0];

nstates = 2;

A=DebrisJacobian(rdeb,adeb,rho0,H,r0,Cd,rhodeb,mu);

phi=reshape(y((nstates+1):end,1),nstates,nstates);
phidot=A*phi;

ydot=NaN(1,nstates+nstates^2);
if adeb>6000
    ydot(1:nstates)=ydot_state;
    ydot((nstates+1):end)=reshape(phidot,1,nstates^2);
end

end