function [J,varargout] = input_fun_SS_analytic_cons_slot_prop(x)

warning('off')

% global VAR
load VAR;

%% extract launch rate

if strcmp(VAR.launch_rate,'generic')
    lam1_S  = x(1:VAR.N_shell);
    lam1_S  = lam1_S(:);
    lam1_Su = x(VAR.N_shell+1:2*VAR.N_shell);
    lam1_Su = lam1_Su(:);
elseif strcmp(VAR.launch_rate,'distribution')
    lam1_S  = x(1)*exp(-(VAR.R02(2:end)-x(2)).^2/x(3)^2)';
    lam1_Su = x(4)*exp(-(VAR.R02(2:end)-x(5)).^2/x(6)^2)';
end

%% propagation

S = sym('S',[VAR.N_shell,1]);
D = sym('D',[VAR.N_shell,1]);
N = sym('N',[VAR.N_shell,1]);
Su = sym('Su',[VAR.N_shell,1]);
lam_S = sym('lam_S',[VAR.N_shell,1]);
lam_Su = sym('lam_Su',[VAR.N_shell,1]);

% compute derivatives

for k=1:VAR.N_shell

    %n=(N(k)+D(k));

    %     F107=(120-70)/2*cos(2*pi*(time/11))+2*40;
    %     F107=(VAR.F107_max-VAR.F107_min)/2*cos(2*pi*(time/11))+VAR.F107_center;
    %
    %     T = 900 + 2.5 *( F107 - 70 )* + 1.5* VAR.Ap;
    %     m = 27 - 0.012* ( VAR.R02(k:k+1) - 200 ); %180 < h(km) < 1000 (is the model valid above 1000 km? if yes, use this, otherwise uncomment the lines below)
    %     if VAR.R02(k)>1000
    %        m(1) = 27 - 0.012* ( 1000 - 200 );
    %     end
    %     if VAR.R02(k+1)>1000
    %        m(2) = 27 - 0.012* ( 1000 - 200 );
    %     end
    %     H = T ./ m;

    % T = 900 + 2.5 *( VAR.F107_middle - 70 ) + 1.5* VAR.Ap_middle;
    % m = 27 - 0.012* ( VAR.R02(k:k+1) - 200 ); %180 < h(km) < 1000 (is the model valid above 1000 km? if yes, use this, otherwise uncomment the lines below)
    % H = T ./ m;

    if k<VAR.N_shell

        n0=(N(k+1));
        n_upper=D(k+1);

        %             rho_k1= 6e-10* exp ( - ( VAR.R02(k+1) - 175 ) / H(2) ); % to include solar flux
        rho_k1 = densityexp(VAR.R02(k+1)); % no solar flux
        rvel_upper_D=-rho_k1*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);%drag flux

    else

        n_upper = 0;
        n0 = 0;
%         m = 27 - 0.012* ( (VAR.R02(k+1)+(VAR.R02(k+1)-VAR.R02(k))) - 200 ); %180 < h(km) < 1000 (is the model valid above 1000 km? if yes, use this, otherwise uncomment the lines below)
%         H_top = T ./ m;
        %             rho_k1= 6e-10* exp ( - ( (VAR.R02(k+1)+(VAR.R02(k+1)-VAR.R02(k))) - 175 ) / H_top ); % to include solar flux
        rho_k1 = densityexp((VAR.R02(k+1)+(VAR.R02(k+1)-VAR.R02(k)))); % no solar flux
        rvel_upper_D=-rho_k1*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k+1)+VAR.R0(k+1)-VAR.R0(k)))*(24*3600*365.25);% drag flux
        rvel_upper_N=-rho_k1*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k+1)+VAR.R0(k+1)-VAR.R0(k)))*(24*3600*365.25);%drag flux

    end

    %         rhok= (6e-10)* exp ( - ( VAR.R02(k) - 175 ) / H(1) ) ; % to include solar flux
    rhok = densityexp(VAR.R02(k)); % no solar flux
    rvel_current_D=-rhok*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);
    rvel_current_N=-rhok*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k)))*(24*3600*365.25);

    % use different probability of collision and generated number of
    % fragments according to the species colliding

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

    A0 = lam_S(k);
    A1 = -1/VAR.Dt;
    A2 = -(VAR.delta+VAR.alpha)*phi_SN;
    A3 = -(VAR.delta+VAR.alpha)*phi_SD;
    A4 = -(1-VAR.slot_par)*VAR.alpha_active*phi_SS;
    A5 = -VAR.alpha_active*phi_SSu*S(k)/(Su(k)+S(k)); 

    Eq1(k) = A0+A1*S(k)+A2*S(k)*N(k)+A3*S(k)*D(k)+A4*S(k)^2 +A5*S(k)*Su(k);
%     Eq1(k) = 0;

    B0 = +n_upper*rvel_upper_D/VAR.Dhu;
    B1 = (1-VAR.P)/VAR.Dt;%S(k)
    B2 = rvel_current_D/VAR.Dhl;%D(k)
    B3 = +VAR.delta*phi_SD;%D(k)*S(k)
    B4 = -phi_DN;%(N(k))*D(k)
    B5 = -phi_DD;%(D(k))*D(k)
    B6 = +VAR.delta*phi_SN;%N(k)*S(k)
    B7 = (1-VAR.P)/VAR.Dt;%Su(k)
    B8 = +VAR.delta*phi_DSu;%D(k)*Su(k)
    B9 = +VAR.delta*phi_NSu;%N(k)*Su(k)

    Eq2(k) = B1*S(k)+B3*D(k)*S(k)+B6*S(k)*N(k)+B4*D(k)*N(k)+B5*D(k)^2+B0+B2*D(k) +B7*Su(k)+B8*D(k)*Su(k)+B9*N(k)*Su(k);

    C0 = +n0*rvel_upper_N/VAR.Dhu;
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

    Eq3(k) = C5*N(k)*S(k)+C2*D(k)*S(k)+C3*N(k)*D(k)+C4*D(k)^2+C6*S(k)^2+C7*N(k)^2+C0+C1*N(k) +C8*D(k)*Su(k)+C9*N(k)*Su(k)+C10*Su(k)^2+C11*S(k)*Su(k);

    A0_Su = lam_Su(k);
    A1_Su = -1/VAR.Dt;
    A2_Su = -(VAR.delta+VAR.alpha)*phi_NSu;
    A3_Su = -(VAR.delta+VAR.alpha)*phi_DSu;
    A4_Su = -VAR.alpha_active*phi_SuSu;
    A5_Su = -VAR.alpha_active*phi_SSu*Su(k)/(Su(k)+S(k)); 

    Eq4(k) = A0_Su+A1_Su*Su(k)+A2_Su*Su(k)*N(k)+A3_Su*Su(k)*D(k)+A4_Su*Su(k)^2 +A5_Su*S(k)*Su(k);

end

var=[S;D;N;Su];
f4=[Eq1.';Eq2.';Eq3.';Eq4.'];
fun4= matlabFunction(f4,'Vars',{var,lam_S,lam_Su});
func_temp = @(x) (fun4(x,lam1_S,lam1_Su));
func=@(t,x) func_temp(x);
options_ode = odeset('reltol', 1.e-8,'abstol', 1.e-8);
[T,X] = ode45(func, VAR.tspan,VAR.x0, options_ode);

S = X(:,1:VAR.N_shell);
D = X(:,VAR.N_shell+1:2*VAR.N_shell);
N = X(:,2*VAR.N_shell+1:3*VAR.N_shell);
Su = X(:,3*VAR.N_shell+1:4*VAR.N_shell);
% N_tot = S+D+N+Su; % total population at each time instant for each shell
% dens_shell = N_tot./(VAR.V*1e-9); % total spatial density over time along all the shells [obj/km^3]
% dens_tot_in = N_tot(1,:)./(VAR.V*1e-9); % Initial total spatial density along all the shells [obj/km^3]
% dens_tot_end = N_tot(end,:)./(VAR.V*1e-9); % Final total spatial density along all the shells [obj/km^3]
N_tot_real_S = sum(S(end,:)); % total number of slotted satellites 
N_tot_real_Su = sum(Su(end,:)); % total number of unslotted satellites 
N_tot_real = N_tot_real_S+N_tot_real_Su; % total number of satellites 
N_diff_debris = (sum(N(end,:))-sum(N(1,:)))/sum(N(1,:)); % difference in the number of debris
% dens_diff_debris = (N(end,:)-N(1,:))./(VAR.V*1e-9); % difference in the number of debris [N/km^3]

N_der = diff(N(end-1:end,:))/(T(end)-T(end-1));

%% propagation time constraint (if the final desired time is not reached)

if T(end)<VAR.tf % propagation doesn't reach the desired final time
    pen_time = 1e9*(VAR.tf-T(end));
else
    pen_time = 0;
end

%% intrinsic capacity constraint

diff_S = VAR.N_sat-S; % if <0, the number of satellites exceeds the maximum allowed
diff_Su = VAR.N_sat-Su; % if <0, the number of satellites exceeds the maximum allowed
pen_sat_S = 0;
pen_const_S = 0;
pen_sat_Su = 0;
pen_const_Su = 0;
% pen_sat_ok = 0;
for i=1:size(X,1) % check for intrinsic capacity constraint for each altitude and for each times step
    pen_sat_it_neg_S = diff_S(i,diff_S(i,:)<0);
%     pen_sat_it_neg_Su = diff_Su(i,diff_Su(i,:)<0);
%     pen_sat_it_pos = diff_S(i,diff_S(i,:)>0);
    if ~isempty(pen_sat_it_neg_S)
        pen_sat_S = pen_sat_S + sum(abs(pen_sat_it_neg_S)); % it penalizes when intrinsic capacity is exceeded
        pen_const_S = pen_const_S + length(VAR.tspan)/i; % it penalizes the time step when intrinsic capacity is exceeded
    end
%     if ~isempty(pen_sat_it_neg_Su)
%         pen_sat_Su = pen_sat_Su + sum(abs(pen_sat_it_neg_Su)); % it penalizes when intrinsic capacity is exceeded
%         pen_const_Su = pen_const_Su + length(VAR.tspan)/i; % it penalizes the time step when intrinsic capacity is exceeded
%     end
%     pen_sat_ok = pen_sat_ok + sum(abs(pen_sat_it_pos)); % NOT USED FOR NOW   
end

% diff_S = VAR.N_sat-S-Su; % if <0, the number of satellites exceeds the maximum allowed
% pen_sat_S = 0;
% pen_const_S = 0;
% pen_sat_Su = 0;
% pen_const_Su = 0;
% % pen_sat_ok = 0;
% for i=1:size(X,1) % check for intrinsic capacity constraint for each altitude and for each times step
%     pen_sat_it_neg_S = diff_S(i,diff_S(i,:)<0);
% %     pen_sat_it_pos = diff_S(i,diff_S(i,:)>0);
%     if ~isempty(pen_sat_it_neg_S)
%         pen_sat_S = pen_sat_S + sum(abs(pen_sat_it_neg_S)); % it penalizes when intrinsic capacity is exceeded
%     end
% %     pen_sat_ok = pen_sat_ok + sum(abs(pen_sat_it_pos)); % NOT USED FOR NOW   
% end


%% final debris population constraint

% if N_diff_debris>VAR.tol_debris % propagation doesn't reach the desired final time
%     pen_debris = 100*N_diff_debris;
% else
%     pen_debris = 0;
% end

N_der_pos = N_der>VAR.tol_debris_der;
if isempty(N_der_pos)
    pen_debris = 0;
else
    pen_debris = 100*sum((N_der(N_der_pos))); % 100
end

%% slotted vs unslotted ratio constraint

% check what happens if this constraint is not included

% if N_tot_real_Su>VAR.threshold*N_tot_real_S 
%     pen_unsl = VAR.beta3*(N_tot_real_Su-VAR.threshold*N_tot_real_S);
% else
%     pen_unsl = 0;
% end

check_S_Su = Su(end,:)-VAR.threshold*S(end,:);
check_S_Su_excess = check_S_Su(check_S_Su<0); %>

if ~isempty(check_S_Su_excess) 
    pen_unsl = VAR.beta3*sum(abs(check_S_Su_excess));
else
    pen_unsl = 0;
end

%% cost function

% max_sat = N_tot_real; % + VAR.beta3*(N_tot_real_Su-VAR.threshold*N_tot_real_S);

J = VAR.beta1/N_tot_real+pen_sat_S+pen_const_S+pen_time+pen_debris+...
    pen_unsl+pen_sat_Su+pen_const_Su;

if nargout==4
    varargout{1} = T;
    varargout{2} = X;
    varargout{3} = [pen_time pen_unsl pen_sat_S+pen_const_S+pen_sat_Su+pen_const_Su pen_debris];
end


end

