% MOCAT4S.m - Model for Collision Assessment Tool with 4 Species
% Models evolution of slotted satellites, unslotted satellites, debris, derelicts
function [OUT] = MOCAT4S(tspan, x0, lam, VAR)
% Inputs:
%   - tspan = times for integration (in years)
%   - x0 = current population VAR.N_shell x 4 in order S, D, N,  Su
%   - lam matrix S, Su per bin so VAR.N_shell x 2 in order S Su
%   - VAR - struct with parameters like mu, beta, phi, etc.
% Outputs:
%   - OUT - struct with propagated state and some other useful info
warning('off')

% disp('Propagating state...')

%% propagation
% Define symbolic variables for populations in each shell k
S = sym('S',[VAR.N_shell,1]);
D = sym('D',[VAR.N_shell,1]);
N = sym('N',[VAR.N_shell,1]);
Su = sym('Su',[VAR.N_shell,1]);
% Symbolic variables for launch rates
lam_S = sym('lam_S',[VAR.N_shell,1]);
lam_Su = sym('lam_Su',[VAR.N_shell,1]);
% Extract numerical values
lam1_S = lam(:,1);
lam1_Su = lam(:,2);
% Initialize accumulator for disposed derelict objects
disposed_derelicts = 0;

% compute derivatives
% Loop over shells
% for k=1:VAR.N_shell % This was the original MOCAT4S ordering.
for k=VAR.N_shell:-1:1 % New ordering to accumulate disposed satellites from higher shells on down
    % Flux terms due to drag
    if k<VAR.N_shell
        % Flux from shell k+1
        n0=(N(k+1));
        n_upper=D(k+1);

        rho_k1 = densityexp(VAR.R02(k+1)); % no solar flux
        rvel_upper_D=-rho_k1*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);%drag flux
    else
        % ASSUMPTION: No flux coming down from highest shell.
        n_upper = 0; 
        n0 = 0;
        rho_k1 = densityexp((VAR.R02(k+1)+(VAR.R02(k+1)-VAR.R02(k)))); % no solar flux
        rvel_upper_D=-rho_k1*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k+1)+VAR.R0(k+1)-VAR.R0(k)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k+1)+VAR.R0(k+1)-VAR.R0(k)))*(24*3600*365.25);%drag flux

    end
    % Flux within shell k
    rhok = densityexp(VAR.R02(k)); % no solar flux
    rvel_current_D=-rhok*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);
    rvel_current_N=-rhok*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);

    % use different probability of collision and generated number of
    % fragments according to the species colliding
    %% Calculate collision probabilities phi_ij
    %% phi_ij is probability of collision between species i and j
    %% Calculated using kinetic theory as phi_ij = pi * v_r * (r_i + r_j)^2 / V(k)
    %% where v_r is relative velocity, r_i and r_j are radii, V(k) is shell volume
    phi_SS = VAR.phi(1,1)/VAR.V(k);
    phi_SD = VAR.phi(1,2)/VAR.V(k);
    phi_SN = VAR.phi(1,3)/VAR.V(k);
    phi_SSu = VAR.phi(1,4)/VAR.V(k);
    phi_DD = VAR.phi(2,2)/VAR.V(k);
    phi_DN = VAR.phi(2,3)/VAR.V(k);
    phi_DSu = VAR.phi(2,4)/VAR.V(k);
    phi_NN = VAR.phi(3,3)/VAR.V(k);
    phi_NSu = VAR.phi(3,4)/VAR.V(k);
    phi_SuSu = VAR.phi(4,4)/VAR.V(k);
    %% Fragmentation coefficients K0_ij
    %% K0_ij is number of fragments generated from collision of i-j
    %% Calculated from NASA breakup model
    K0_SS = VAR.K0(1,1);
    K0_SD = VAR.K0(1,2);
    K0_SN = VAR.K0(1,3);
    K0_SSu = VAR.K0(1,4);
    K0_DD = VAR.K0(2,2);
    K0_DN = VAR.K0(2,3);
    K0_DSu = VAR.K0(2,4);
    K0_NN = VAR.K0(3,3);
    K0_NSu = VAR.K0(3,4);
    K0_SuSu = VAR.K0(4,4);
    % Derivatives for species in shell k
    A0 = lam_S(k); % A0 - launch rate
    A1 = -1/VAR.Dt; % A1 - decay rate
    A2 = -(VAR.delta+VAR.alpha)*phi_SN; % A2 - collisions with debris
    A3 = -(VAR.delta+VAR.alpha)*phi_SD; % A3 - collisions with derelicts
    A4 = -(1-VAR.slot_par)*VAR.alpha_active*phi_SS; % A4 - collisions between slotted satellites
    A5 = -VAR.alpha_active*phi_SSu*S(k)/(Su(k)+S(k)); % A5 - collisions with unslotted satellites

    % Calculate time derivative for slotted satellites in shell k
    Eq1(k) = A0 + A1*S(k) + A2*S(k)*N(k) + A3*S(k)*D(k) + A4*S(k)^2 + A5*S(k)*Su(k);
    % Terms for:
    % - Launch
    % - Decay
    % - Collisions with other species
    % - Flux from higher shell

    B0 = +n_upper*rvel_upper_D/VAR.Dhu; % flux from above due to drag
    B1 = (1-VAR.P)/VAR.Dt;%S(k) % creation from retired slotted satellites
    B2 = rvel_current_D/VAR.Dhl;%D(k) % decay due to drag
    B3 = +VAR.delta*phi_SD;%D(k)*S(k) % collisions S-D
    B4 = -phi_DN;%(N(k))*D(k) % collisions D-N 
    B5 = -phi_DD;%(D(k))*D(k) % collisions D-D
    B6 = +VAR.delta*phi_SN;%N(k)*S(k) % collisions N-S
    B7 = (1-VAR.P)/VAR.Dt;%Su(k) % creation from retired unslotted satellites
    B8 = +VAR.delta*phi_DSu;%D(k)*Su(k) % collisions D-Su 
    B9 = +VAR.delta*phi_NSu;%N(k)*Su(k) % collisions N-Su

    % Calculate time derivative for derelicts in shell k
    %% If the current index is within the naturally-compliant region. Default MOCAT4S behavior, no change.
    if k < VAR.k_star
        Eq2(k) = B1*S(k) + B3*D(k)*S(k) + B6*S(k)*N(k) + B4*D(k)*N(k) + B5*D(k)^2 + B0 + B2*D(k) + B7*Su(k) + B8*D(k)*Su(k) + B9*N(k)*Su(k);
    end
    % Terms for:
    % - Creation from slotted satellite retirement 
    % - Collisions with other species
    % - Atmospheric drag
    % - Flux from higher shell
    %% If the current index is not within the naturally-compliant region. Take accumulation from B1 and B7 and store it in disposed_derelicts
    if k > VAR.k_star
        Eq2(k) = B3*D(k)*S(k) + B6*S(k)*N(k) + B4*D(k)*N(k) + B5*D(k)^2 + B0 + B2*D(k) + B8*D(k)*Su(k) + B9*N(k)*Su(k);
        disposed_derelicts = disposed_derelicts + B1*S(k) + B7*Su(k);
    end
    %% If the current shell is the boundary of the naturally-compliant region. Inject disposed_derelicts into here
    if k == VAR.k_star
        Eq2(k) = B1*S(k) + B3*D(k)*S(k) + B6*S(k)*N(k) + B4*D(k)*N(k) + B5*D(k)^2 + B0 + B2*D(k) + B7*Su(k) + B8*D(k)*Su(k) + B9*N(k)*Su(k) + disposed_derelicts; 
    end

    C0 = +n0*rvel_upper_N/VAR.Dhu; % flux from above due to drag
    C1 = rvel_current_N/VAR.Dhl;%N(k)
    C2 = +phi_SD*VAR.alpha*K0_SD;%D(k)*S(k) %new
    C3 = K0_DN*phi_DN;%(N(k))*D(k)
    C4 = K0_DD*phi_DD;%(D(k))*D(k)
    C5 = phi_SN*VAR.alpha*K0_SN;%(N(k))*S(k)
    C6 = (1-VAR.slot_par)*VAR.alpha_active*phi_SS*K0_SS;%S(k)^2
    C7 = phi_NN*K0_NN;%N(k)^2 %new
    C8 = +phi_DSu*VAR.alpha*K0_DSu;%D(k)*Su(k) %new
    C9 = phi_NSu*VAR.alpha*K0_NSu;%(N(k))*S(k)
    C10 = VAR.alpha_active*phi_SuSu*K0_SuSu;
    C11 = +VAR.alpha_active*K0_SSu*phi_SSu; % S(k)*Su(k) 

    % Calculate time derivative for debris in shell k 
    Eq3(k) = C5*N(k)*S(k) + C2*D(k)*S(k) + C3*N(k)*D(k) + C4*D(k)^2 + C6*S(k)^2 + C7*N(k)^2+C0 + C1*N(k) + C8*D(k)*Su(k) + C9*N(k)*Su(k) + C10*Su(k)^2 + C11*S(k)*Su(k);
    % Terms for:
    % - Creation from collisions 
    % - Atmospheric drag
    % - Flux from higher shell

    A0_Su = lam_Su(k); % Launch rate of unslotted satellites
    A1_Su = -1/VAR.Dt; % Decay rate of unslotted satellites
    A2_Su = -(VAR.delta+VAR.alpha)*phi_NSu; % Collisions with debris
    A3_Su = -(VAR.delta+VAR.alpha)*phi_DSu; % Collisions with derelicts
    A4_Su = -VAR.alpha_active*phi_SuSu; % Collisions between unslotted satellites
    A5_Su = -VAR.alpha_active*phi_SSu*Su(k)/(Su(k)+S(k)); % Collisions with slotted satellites

    % Calculate time derivative for unslotted satellites in shell k
    Eq4(k) = A0_Su + A1_Su*Su(k) + A2_Su*Su(k)*N(k) + A3_Su*Su(k)*D(k) + A4_Su*Su(k)^2 + A5_Su*S(k)*Su(k);
    % Terms for:
    % - Launch
    % - Decay
    % - Collisions with other species
    % - Flux from higher shell

end

% Collect derivative expressions into vector
var=[S;D;N;Su];
f4=[Eq1.';Eq2.';Eq3.';Eq4.'];
% Convert to MATLAB function
fun4= matlabFunction(f4,'Vars',{var,lam_S,lam_Su});
% Function calling with numerical lam
func_temp = @(x) (fun4(x,lam1_S,lam1_Su));
func=@(t,x) func_temp(x);
% Integrate ODE system
options_ode = odeset('reltol', 1.e-8,'abstol', 1.e-8);
[T,X] = ode45(func, tspan, x0, options_ode);

%%% Extract populations
S = X(:,1:VAR.N_shell);
D = X(:,VAR.N_shell+1:2*VAR.N_shell);
N = X(:,2*VAR.N_shell+1:3*VAR.N_shell);
Su = X(:,3*VAR.N_shell+1:4*VAR.N_shell);

% N_tot = S+D+N+Su; % total population at each time instant for each shell
% dens_shell = N_tot./(VAR.V*1e-9); % total spatial density over time along all the shells [obj/km^3]
% dens_tot_in = N_tot(1,:)./(VAR.V*1e-9); % Initial total spatial density along all the shells [obj/km^3]
% dens_tot_end = N_tot(end,:)./(VAR.V*1e-9); % Final total spatial density along all the shells [obj/km^3]
% N_tot_real_S = sum(S(end,:)); % total number of slotted satellites 
% N_tot_real_Su = sum(Su(end,:)); % total number of unslotted satellites 
% N_tot_real = N_tot_real_S+N_tot_real_Su; % total number of satellites 
% N_diff_debris = (sum(N(end,:))-sum(N(1,:)))/sum(N(1,:)); % difference in the number of debris
% dens_diff_debris = (N(end,:)-N(1,:))./(VAR.V*1e-9); % difference in the number of debris [N/km^3]
% N_der = diff(N(end-1:end,:))/(T(end)-T(end-1));

OUT = struct;
OUT.S = S; % slotted satellites
OUT.D = D; % derelicts
OUT.N = N; % debris
OUT.Su = Su; % unslotted satellites
OUT.T = T;

% Print message that state has been propagated
% fprintf('State propagated.\n');

end
