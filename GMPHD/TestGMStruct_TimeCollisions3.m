clear all
close all
clc

% Mark Moretto and Akhil Rao
% September 2023
% Test GM implementation

%% Set time horizon and integration steps
nyears=1; % Number of years to propagate for horizon
ntime=nyears*10; % Number of integration time steps within horizon
dicoll=10; % Number of timesteps between collisions

%% Test -- replace this with loading in the TLE initial conditions
constellation_levels = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
                        0, 0, 0, 0, 0, 0, 0, 0, 0, 0].^0 * 0;
fringe_levels = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0,...
                    0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
constellation_launches = constellation_levels;
constellation_launches(10) = 300;
fringe_launches = fringe_levels;
fringe_launches(8) = 100;

% Load parameter values
GMPHD_params = GMPHD_VAR_Cons();
GMPHD_params.r0 = 700+GMPHD_params.rE; % reference height for drag

% Build shellstruct to hold discretized objects for reporting and collision calculations
shellstruct = buildShellStruct(constellation_levels, fringe_levels, constellation_launches, fringe_launches, GMPHD_params);

initial_GM_objects = prepareInitialGMs(GMPHD_params, nyears, ntime,dicoll);
t = initial_GM_objects.t;
GM = initial_GM_objects.GM;
i_collisions = initial_GM_objects.i_collisions;
GM_coll = initial_GM_objects.GM_coll;

tic
[t,GM1,shellstruct] = PropagateGM_timehist_collisions(t, GM, GMPHD_params, shellstruct, i_collisions, GM_coll);
toc

%% Testing iam integration
for i = 1:40 % Currently hard-coded to match MOCAT shell discretization for launch and reporting; underlying GMs are defined over continuous altitude range. Debris objects not initialized here since they are also defined over continuous radius range.
    Si(i) = shellstruct(i).nconst(end);
    Sui(i) = shellstruct(i).nfringe(end);
    lam(i,1) = shellstruct(i).launch_rate_const(end);
    lam(i,2) = shellstruct(i).launch_rate_fringe(end);
end

figure(2)
plot(t(i_collisions)./(365.25*24*60*60),shellstruct(1).nconst,'c')
hold on
plot(t(i_collisions)./(365.25*24*60*60),shellstruct(1).nfringe,'c.-')

plot(t(i_collisions)./(365.25*24*60*60),shellstruct(2).nconst,'m')
plot(t(i_collisions)./(365.25*24*60*60),shellstruct(2).nfringe,'m.-')
legend('Lower Shell, Const','Lower Shell, Fringe','Upper Shell, Const','Upper Shell, Fringe')
xlabel('Time [Years]')
ylabel('Num Sats')
% save figure to file named GM2.png
print -dpng GM2.png

% Plot collisions at each time step in each shell
figure(3)
for i=1:16
    plot(t(i_collisions)./(365.25*24*60*60),shellstruct(i).ncoll)
    hold on
end
xlabel('Time [Years]')
ylabel('Number Collisions')
legend('Shell 1','Shell 2','Shell 3','Shell 4','Shell 5','Shell 6','Shell 7','Shell 8','Shell 9','Shell 10','Shell 11','Shell 12','Shell 13','Shell 14','Shell 15','Shell 16')
% save figure to file named GM3.png
print -dpng GM3.png

% Post-processing block to calculate debris numbers in discrete categories for reporting
rbounds=[0, 0.1/1000, 0.5/1000, 1000]; % Size categories for debris in post-processing. 0-10cm, 10-50cm, 50cm-1m, 1m+
abounds=[200:35:1600]+GMPHD_params.rE; % altitude bins for debris in post-processing. 300-1600km

%% Numerical integration over Gaussians. Runs in parallel.
tic
[D_lnt, D_st, D_lt, Si, Sui, lam]=PostprocessDebris(GM1, t, abounds, rbounds, GMPHD_params, shellstruct, 1, length(t));
toc

% Extract numdeb values to three vectors: D_lnt, D_st, D_lt. Each is a k slice from (end,:,k) of numdeb
D_lnt = numdeb(end,:,1);
D_st = numdeb(end,:,2);
D_lt = numdeb(end,:,3);
debris_t = [D_lnt', D_st', D_lt'];

for i=1:(length(abounds)-1)
figure(3+i)
plot(t./(365.25*24*60*60),squeeze(numdeb(:,i,:)))
xlabel("Time [years]")
ylabel("Num Debris []")
legend("Lethal Non-Trackable","Trackable","Large")
title(strjoin(["Shell",num2str(abounds(i)-GMPHD_params.rE),"-",num2str(abounds(i+1)-GMPHD_params.rE),"km Altitude"]))
end



figure(2)
plot(t(i_collisions)./(365.25*24*60*60),shellstruct(1).nconst,'c')
hold on
plot(t(i_collisions)./(365.25*24*60*60),shellstruct(1).nfringe,'c.-')

plot(t(i_collisions)./(365.25*24*60*60),shellstruct(2).nconst,'m')
plot(t(i_collisions)./(365.25*24*60*60),shellstruct(2).nfringe,'m.-')
legend('Lower Shell, Const','Lower Shell, Fringe','Upper Shell, Const','Upper Shell, Fringe')
xlabel('Time [Years]')
ylabel('Num Sats')
% save figure to file named GM2.png
print -dpng GM2.png

% Plot collisions at each time step in each shell
figure(3)
for i=1:16
    plot(t(i_collisions)./(365.25*24*60*60),shellstruct(i).ncoll)
    hold on
end
xlabel('Time [Years]')
ylabel('Number Collisions')
legend('Shell 1','Shell 2','Shell 3','Shell 4','Shell 5','Shell 6','Shell 7','Shell 8','Shell 9','Shell 10','Shell 11','Shell 12','Shell 13','Shell 14','Shell 15','Shell 16')
% save figure to file named GM3.png
print -dpng GM3.png


rbounds=[0,0.1/1000,0.5/1000,3/1000]
abounds=[300:35:1700]+GMPHD_params.rE;

[numdeb]=PostprocessDebris(GM1,t,abounds,rbounds);

for i=1:(length(abounds)-1)
figure(100+i)
plot(t./(365.25*24*60*60),squeeze(numdeb(:,i,:)))
xlabel("Time [years]")
ylabel("Num Debris []")
legend("Lethal Non-Trackable","Trackable","Large")
title(strjoin(["Shell",num2str(i),num2str(abounds(i)-GMPHD_params.rE),"-",num2str(abounds(i+1)-GMPHD_params.rE),"km Altitude"]))
print(strjoin(["Shell",num2str(i),num2str(abounds(i)-GMPHD_params.rE),"-",num2str(abounds(i+1)-GMPHD_params.rE),"km_Altitude.png"],''),'-dpng')
end

%%

nx=1000
ny=100

xim=linspace(GMPHD_params.rE,GMPHD_params.rE+1500,nx);
yim=linspace(0,3/1000,ny);

v = VideoWriter('debris_10year.mp4','MPEG-4');
open(v)
for i = 1:length(t)
    figure(1000)
    [h,X1,X2,yi0] = PlotDist(xim,yim,squeeze(GM1(i,:)),1);
    clim([0,3000000])
    colorbar
    title(t(i)/(365.25*24*60*60))
    drawnow
    frame = getframe(gcf);

    writeVideo(v,frame);

end

close(v)

