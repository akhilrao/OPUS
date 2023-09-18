function [cdot]=calcCdGaussian(ndeb,xbar,P,rmin,rmax,rdebmin,rdebmax,V,rsat,nsat,ffail)
% This function computes the collision probability assuming the debris is distributed according to a Gaussian distribution.
% KMS units

    detP=det(P); % determinant of the covariance matrix

    Vshell=4/3*pi*(rmax^3-rmin^3); % volume of the shell

    fint=@(r,rdeb) cdGauss(r,rdeb,ffail,nsat,ndeb,V,rsat,xbar,P,Vshell);

    cdot=integral2(fint,rmin,rmax,rdebmin,rdebmax,'AbsTol',1e-12,'RelTol',1e-12);

end


function [cdint]=cdGauss(r,rdeb,ffail,nsat,ndeb,V,rsat,xbar,P,Vshell)
% This function calculates the collision probability between a satellite and a debris object. Assumes satellites never collide with each other.

    cdint=ffail*nsat*ndeb*V*(rsat^2+rdeb(:).^2)*pi.*mvnpdf([r(:),rdeb(:)],xbar,P)/Vshell;
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