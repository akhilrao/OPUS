function [cdot]=calcNumGaussian(ndeb,xbar,P,rmin,rmax,rdebmin,rdebmax)
% This function computes the collision probability assuming the debris is distributed according to a Gaussian distribution.
% KMS units

    fint=@(r,rdeb) NumGauss(r,rdeb,xbar,P,ndeb);

    cdot=integral2(fint,rmin,rmax,rdebmin,rdebmax,'AbsTol',1e-12,'RelTol',1e-12);

end

function [cdint]=NumGauss(r,rdeb,xbar,P,ndeb)
% This function calculates the collision probability between a satellite and a debris object. Assumes satellites never collide with each other.

    cdint=ndeb*mvnpdf([r(:),rdeb(:)],xbar,P);
    % Elements of the formula:
    % ffail: failure rate
    % nsat: number of satellites
    % ndeb: number of debris
    % V: relative velocity
    % rsat: satellite radius
    % rdeb: debris radius
    % pi: pi
    % mvnpdf: multivariate normal probability density function
    % Vshell: volume of the shell

    cdint=reshape(cdint,size(r));

end