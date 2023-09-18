function [D_lnt, D_st, D_lt, Si, Sui, lam]=PostprocessDebris(GM, t, aarr, rarr, GMPHD_params, shellstruct, pb, n_workers)
% Function to calculate number of debris objects across discrete categories given by user, as well as other post-processing of object stocks and launch rates.
% Inputs:
%   GM: Gaussian mixture model of debris objects
%   t: time vector
%   aarr: array of altitude bins
%   rarr: array of radius bins
%   GMPHD_params: structure containing GMPHD parameters
%   shellstruct: structure containing shell parameters
%   pb: boolean to display progress indicator
%   n_workers: number of workers to use in parallel pool
% Outputs:
%   numdeb: number of debris objects across discrete (time, altitude, radius) categories

% Initialize number of debris tensor
numdeb=NaN(length(t),length(aarr)-1,length(rarr)-1);

if n_workers > 1
    parpool('local', n_workers)
end

% Calculate number of debris across discrete categories set by user. i index is time, j index is altitude, k index is radius
for j = 1:(length(aarr)-1)
    if pb
        disp(['Debris post-processing progress: ', num2str(floor(j/(length(aarr)-1)*100)), '%'])  % progress indicator
    end
    for k = 1:(length(rarr)-1)
        parfor i = 1:length(t)
            GMi=squeeze(GM(i,:));
            rmin=aarr(j);
            rmax=aarr(j+1);
            rdebmin=rarr(k);
            rdebmax=rarr(k+1);
            [num, ~]=calcNumGaussianStruct(GMi,rmin,rmax,rdebmin,rdebmax);
            numdeb(i,j,k)=num;
        end
    end
end

p = gcp('nocreate'); % Get the current parallel pool (if it exists)
if ~isempty(p)
    delete(p); % Shut down the parallel pool
end

% Extract numdeb values to three vectors: D_lnt, D_st, D_lt. Each is a k slice from (end,:,k) of numdeb
D_lnt = numdeb(end,:,1);
D_st = numdeb(end,:,2);
D_lt = numdeb(end,:,3);
debris_t = [D_lnt, D_st, D_lt];

% Create Si, Sui, lam objects to hold initial stocks and launch rates. This doesn't really depend on anything else in this code, it's just here to simplify the calls in iam_solver.m
for i = 1:40 % Currently hard-coded to match MOCAT shell discretization for launch and reporting; underlying GMs are defined over continuous altitude range.
    Si(i) = shellstruct(i).nconst(end);
    Sui(i) = shellstruct(i).nfringe(end);
    lam(i,1) = shellstruct(i).launch_rate_const(end);
    lam(i,2) = shellstruct(i).launch_rate_fringe(end);
end

end