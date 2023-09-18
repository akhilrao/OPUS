function [cdot_total, cdot_ind]=calcNumGaussianStruct(GM,rmin,rmax,rdebmin,rdebmax)
    % Function to wrap calcNumGaussian and calculate number of collisions
    % Inputs: 
    %   GM0: Gaussian mixtures of debris state
    %       rmin: minimum radius of shell
    %       rmax: maximum radius of shell
    %       rdebmin: minimum radius of debris
    %       rdebmax: maximum radius of debris
    %       V: volume of shell
    %       r: radius of shell
    %       n: number of satellites in shell
    %       ffail: failure probability
    %   sat_type: 'fringe' or 'const'
    % Outputs:
    %   cdot_total: total number of collisions
    %   cdot_ind: number of collisions for each Gaussian
cdot_ind=NaN(length(GM),1);

cdot_total=0;

for i=1:length(GM)
    
    if ~isempty(GM(i).state)
    if isfinite(GM(i).state(1))

    [cdot]=calcNumGaussian(GM(i).c,GM(i).state,GM(i).P,rmin,rmax,rdebmin,rdebmax);

    cdot_ind(i,1)=cdot;
    cdot_total=cdot_total+cdot;

    end
    end

end


end