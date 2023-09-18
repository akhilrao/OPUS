function [VAR] = MOCAT4S_VAR_Cons()

% Constants
rad = pi/180;
mu_km=398600.4415;%km^3/s^2
mu=3.986004418e14;%meters^3/s^2
re = 6378.1366; % [km]
years_s = 365*24*3600;

% Parameters needed for all MOCAT-SSEMs
VAR.N_shell = 40;
VAR.h_max = 1600;
VAR.h_min = 200;
VAR.N_step = 5; % how many steps to calculate between t0 and tf
N_shell = VAR.N_shell;
h_max = VAR.h_max;
h_min = VAR.h_min;
R0=linspace(h_min,h_max,N_shell+1);
R02=linspace(h_min,h_max,N_shell+1);
VAR.deltaH=R02(2)-R02(1); % thickness of the shell [km]

% Parameters needed for MOCAT-4S
VAR.alpha = 0.2; % fraction of collisions that satellites fail to avoid.
VAR.alpha_active=0.01; % fraction of collisions occuring between active satellites 
VAR.Dt = 5; % operation lifetime: after those years, active satellites become derelict [years]
VAR.slot_par = 1; % slotting effectiveness parameter in [0,1] where 1 = no cols. between slotted satellites.
VAR.delta = 10; %10
VAR.P = 0; % probability of regulatory non-compliance. in the original MOCAT4S this was "probability of instant disposal" and set to 0.5.
VAR.Cd = 2.2; % Coefficent of drag for drag purposes. unitless
VAR.mass      = [223 223 0.640 223];      % median mass values from MASTER [kg] SOMMA
VAR.A         = [1.741 1.741 0.020 1.741]; % median area values from MASTER [m^2] SOMMA
VAR.diameter  = [1.490 1.490 0.180 1.490]; % median diameter values from MASTER [m] SOMMA
VAR.v_imp = 10; % impact velocity [km/s]
VAR.LC = 0.1; % minimum size of fragments [m]
VAR.deg_intrinsic = 45; 

VAR.sep_dist_method = "distance";% either 'distance' or 'angle'. Sets minimum separation distance for intrinsic capacity constraint.
VAR.sep_angle = 0.1; % minimum angular separation distance [deg]. Only used if sep_dist_method VAR.sep_dist_method = 'angle'
VAR.sep_dist = 5; % minimum separation distance [km]. Only used if sep_dist_method VAR.sep_dist_method = 'distance'

VAR.input_pop = 'TLE';  % 'TLE' or 'distribution'.  TLE pulls binned present day objects. 'distribution' uses an arbitary hardcoded population distribution.
% Generate TLE binned appropriately using the getTLEBins python method
% TODO: 
% savedir = py.tletobins.getTLEBins(VAR.deltaH, VAR.h_min, VAR.h_max, savedir = None, graph = False)
% VAR.filename_N = fullfile(savedir, 'Counts_DEBRIS_bins_' +num2str(VAR.deltaH) + '.csv';
% VAR.filename_S = fullfile(savedir, 'Counts_PAYLOADslot_bins_' +num2str(VAR.deltaH) + '.csv'
% VAR.filename_Su = fullfile(savedir, 'Counts_PAYLOADunslot_bins_' +num2str(VAR.deltaH) + '.csv'
% VAR.filename_D = fullfile(savedir, 'Counts_DERELICT_bins_' +num2str(VAR.deltaH) + '.csv'

% Load population files
VAR.filename_N = 'x0_TLE/Counts_DEBRIS_bins_35.csv';
VAR.filename_S = 'x0_TLE/Counts_PAYLOADslot_bins_35.csv';
VAR.filename_Su = 'x0_TLE/Counts_PAYLOADunslot_bins_35.csv';
VAR.filename_D = 'x0_TLE/Counts_DERELICT_bins_35.csv';

R1=R0;
R0=(re+R0)*1000; % orbital radius in meters
% V=4*pi*R0.^2*deltaH*1000;
VAR.V=4/3*pi*(diff(R0.^3)); % volume of the shells [m^3]
VAR.v=VAR.v_imp*1000*(24*3600*365.25);% impact velocity [m/year]
VAR.Dhl=VAR.deltaH*1000; % Shell thickness parameter [m]
VAR.Dhu=-VAR.deltaH*1000; % Shell thickness parameter [m]
Cd=VAR.Cd;

%% physical characteristics
n_obj = length(VAR.mass);
radii = VAR.diameter/2;
VAR.area_mass = VAR.A./VAR.mass; % area/mass ratio  [# shells, # object types]
beta = VAR.Cd*VAR.area_mass; % ballistic coefficient [# shells, # object types]. 

%% NASA standard breakup models implementation for collisions
v_imp = VAR.v_imp; % impact velocity [km/s]
LC = VAR.LC; % minimum size of fragments [m]
n_f_catastrophic = @(M1,M2) 0.1*LC^(-1.71)*(M1+M2)^(0.75); % number of fragments generated during a catastrophic collision (NASA standard break-up model). M is the sum of the mass of the objects colliding in kg
n_f_damaging = @(M1,M2) 0.1*LC^(-1.71)*(min(M1,M2)*v_imp^2)^(0.75); % number of fragments generated during a non-catastrophic collision (improved NASA standard break-up model: takes into account the kinetic energy). M is the mass of the less massive object colliding in kg
K0 = zeros(n_obj,n_obj); % number of fragments generated in each collisions among the species for each altitude bin
ind_obj_catastrophic = 3; % hypothesis of catastrophic collisions only between satellites and derelict satellites
sigma = zeros(n_obj,n_obj);
for q = 1:length(radii)
    mass_q = VAR.mass(q);
    for qq = 1:length(radii)
        sigma(q,qq) = (radii(q)+radii(qq))^2;
        mass_qq = VAR.mass(qq);
        if q~=ind_obj_catastrophic && qq~=ind_obj_catastrophic
            K0(q,qq) = n_f_catastrophic(mass_q,mass_qq); % catastrophic collision
        else
            K0(q,qq) = n_f_damaging(mass_q,mass_qq); % (non-catastrophic) damaging collision
        end
    end
end

% phi=sigma*v; % adopted before 
f_corrective = 1; %0.265 this is a corrective factor (obtained by try-and-error) that led to obtain a similar number of fragments wrt Somma
% phi=f_corrective*1/2*pi*sigma*v; % new (from Somma) the factor 1/2 is not correct
phi=f_corrective*pi*sigma*VAR.v; % new (from Somma)

%% load results from intrinsic capacity and inclination vs altitude
% Note that this currently not used because there is no constrained
% optimization.
% TODO: Add warning to actual model execution if S exceeds N_sat_tot_intrinsic?
filename = 'MOCAT4S/LawDF_min_500_fit_range_max_10000_sols_1.csv';
A = readmatrix(filename);
A = A(2:end,:);
index_mat = A(:,1);
inc_mat   = A(:,2); % [deg]
c_mat     = A(:,3);
b_mat     = A(:,4);
R2_mat    = A(:,5);
[minDistance, ind_intrinsic]=min(abs(VAR.deg_intrinsic-A(:,2))); %Pull the closest interp from the table
if minDistance > 0.1
    warning("The interpolated intrinsic capacity powerlaw is more than 0.1 degrees from the chosen inclination and may not be accurate. For retrograde orbits, use prograde equivalent.")
end
% sep_dist_method = 'angle'; % choose angle or distance for the minimum safe distance
sep_dist_method = VAR.sep_dist_method;

if strcmp(sep_dist_method,'angle')
    sep_angle = VAR.sep_angle; % minimum angular separation distance [deg]
    N_sat_eq_ang = @(c,b) (sep_angle./c).^(1./b); % results from intrinsic capacity analysis
    N_sat = N_sat_eq_ang(c_mat(ind_intrinsic),b_mat(ind_intrinsic)).*(R02(2)-R02(1))./(sep_angle*rad*(R02(2:end)+re)); % N_sat_intrinsic*(n° shells per bin)
else
    sep_dist = VAR.sep_dist; % minimum separation distance [km], 5
    N_sat_eq = @(h,c,b) (sep_dist./(re+h)./rad./c).^(1./b); % results from intrinsic capacity analysis
    N_sat = N_sat_eq(R02(2:end),c_mat(ind_intrinsic),b_mat(ind_intrinsic))*(R02(2)-R02(1))/sep_dist; % N_sat_intrinsic*(n°  shells per bin)
end

N_sat_tot_intrinsic = sum(N_sat);

%% initial conditions

%input_pop = 'TLE';
% input_pop = 'distribution';

if strcmp(VAR.input_pop,'TLE')
    %%% MODIFY FOR MORE SHELLS
    % load initial conditions from TLEs
    
    filename_N = VAR.filename_N;
    filename_S = VAR.filename_S;
    filename_Su = VAR.filename_Su;
    filename_D = VAR.filename_D;

    Ni = readmatrix(filename_N);
    Ni = Ni(:,2)';
    
    Si = readmatrix(filename_S);
%   Si = Si(:,3)';
    Si = Si(:,2)';
       
    Sui = readmatrix(filename_Su);
    Sui = Sui(:,2)';
       
    Di = readmatrix(filename_D);
    Di = Di(:,2)';
%   Di = floor(0.5*(Si+Sui));

elseif strcmp(VAR.input_pop,'distribution')

    % previous initial conditions

    mu_S=500;sig_S=300^2;S0=1000;
    mu_Su=650;sig_Su=300^2;Su0=850;
    mu_D=450;sig_D=300^2;D0=500;
    mu_N=300;sig_N=150^2;N0=400;
    mu_lam_S=350;sig_lam_S=300^2;lam0_S=1000; %lam0=3000
    mu_lam_Su=350;sig_lam_Su=300^2;lam0_Su=1000; %lam0=3000
    lam_model_S=lam0_S*exp(-(R02(2:end)-mu_lam_S).^2/sig_lam_S); 
    lam_model_Su=lam0_Su*exp(-(R02(2:end)-mu_lam_Su).^2/sig_lam_Su); 

    Si=S0*exp(-(R02(2:end)-mu_S).^2/sig_S);
    Sui=Su0*exp(-(R02(2:end)-mu_Su).^2/sig_Su);
    Di=D0*exp(-(R02(2:end)-mu_D).^2/sig_D);
    Ni=N0*exp(-(R02(2:end)-mu_N).^2/sig_N);

end

%% start computation


options = odeset('reltol', 1.e-4,'abstol', 1.e-4);

VAR.x0 = [Si';Di';Ni';Sui'];
VAR.options = options;
%VAR.h = h;
VAR.beta = beta;
VAR.mu = mu;
VAR.mu_km = mu_km;
VAR.K0 = K0;
VAR.phi = phi;
VAR.R0 = R0;
VAR.R02 = R02;
VAR.N_sat = N_sat;

launch_rate = 'generic';
VAR.S_Su = 0.1; %ratio of unslotted to slotted spacecraft with generic or distribution model
MOCAT = '4S'; 
VAR.launch_rate = launch_rate;

