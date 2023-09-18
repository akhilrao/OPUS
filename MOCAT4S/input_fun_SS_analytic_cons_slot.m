function [J,varargout] = input_fun_SS_analytic_cons_slot(x)

warning('off')

% global VAR
load VAR;

if strcmp(VAR.launch_rate,'generic')
    lam_S  = x(1:VAR.N_shell);
    lam_S  = lam_S(:);
    lam_Su = x(VAR.N_shell+1:2*VAR.N_shell);
    lam_Su = lam_Su(:);
elseif strcmp(VAR.launch_rate,'distribution')
    lam_S  = x(1)*exp(-(VAR.R02(2:end)-x(2)).^2/x(3)^2)';
    lam_Su = x(4)*exp(-(VAR.R02(2:end)-x(5)).^2/x(6)^2)';
end

S_eq = lam_S*VAR.Dt;
Su_eq = lam_Su*VAR.Dt;

k = VAR.k;
n0=VAR.n0;
n_upper=VAR.n_upper;
try

    syms S D N Su
    % compute derivatives

    
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

        T = 900 + 2.5 *( VAR.F107_middle - 70 ) + 1.5* VAR.Ap_middle;
        m = 27 - 0.012* ( VAR.R02(k:k+1) - 200 ); %180 < h(km) < 1000 (is the model valid above 1000 km? if yes, use this, otherwise uncomment the lines below)
        H = T ./ m;

        if k<VAR.N_shell

%             rho_k1= 6e-10* exp ( - ( VAR.R02(k+1) - 175 ) / H(2) ); % to include solar flux
            rho_k1 = densityexp(VAR.R02(k+1)); % no solar flux
            rvel_upper_D=-rho_k1*VAR.beta(2)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);% drag flux
            rvel_upper_N=-rho_k1*VAR.beta(3)*sqrt(VAR.mu*(VAR.R0(k+1)))*(24*3600*365.25);%drag flux

        else

            n_upper = D;
            n0 = N;            
            m = 27 - 0.012* ( (VAR.R02(k+1)+(VAR.R02(k+1)-VAR.R02(k))) - 200 ); %180 < h(km) < 1000 (is the model valid above 1000 km? if yes, use this, otherwise uncomment the lines below)
            H_top = T ./ m;
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

        A0 = lam_S;
        A1 = -1/VAR.Dt;
        A2 = -(VAR.delta+VAR.alpha)*phi_SN;
        A3 = -(VAR.delta+VAR.alpha)*phi_SD;
        A4 = -(1-VAR.slot_par)*VAR.alpha_active*phi_SS;
        A5 = -alpha_active*phi_SSu*S/(Su+S);

        Eq1 = A0+A1*S+A2*S*N+A3*S*D+A4*S^2 +A5*S*Su;

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

        Eq2 = B1*S+B3*D*S+B6*S*N+B4*D*N+B5*D^2+B0+B2*D +B7*Su+B8*D*Su+B9*N*Su;

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
        C10 = (1-VAR.slot_par)*VAR.alpha_active*phi_SuSu*K0_SuSu;%Su(k)^2
        C11 = +alpha_active*K0_SSu*phi_SSu; % S(k)*Su(k)

        Eq3 = C5*N*S+C2*D*S+C3*N*D+C4*D^2+C6*S^2+C7*N^2+C0+C1*N +C8*D*Su+C9*N*Su+C10*Su^2+C11*S*Su;

        A0_Su = lam_Su;
        A1_Su = -1/VAR.Dt;
        A2_Su = -(VAR.delta+VAR.alpha)*phi_NSu;
        A3_Su = -(VAR.delta+VAR.alpha)*phi_DSu;
        A4_Su = -(1-VAR.slot_par)*VAR.alpha_active*phi_SuSu;
        A5_Su = -alpha_active*phi_SSu*Su/(Su+S);

        Eq4 = A0_Su+A1_Su*Su+A2_Su*Su*N+A3_Su*Su*D+A4_Su*Su^2 +A5_Su*S*Su;

        out=solve([Eq1;Eq2;Eq3;Eq4]);

        sol_temp = double([vpa(out.S) vpa(out.D) vpa(out.N) vpa(out.Su)]);
        sol_temp = sol_temp((abs(imag(sol_temp(:,1)))<VAR.thres_imag & ...
            abs(imag(sol_temp(:,2)))<VAR.thres_imag & abs(imag(sol_temp(:,3)))<VAR.thres_imag & abs(imag(sol_temp(:,4)))<VAR.thres_imag),:);
        [~,ind_m] = max(sol_temp(sol_temp(:,1)>0 & sol_temp(:,2)>0 & sol_temp(:,3)>0 & sol_temp(:,4)>0,1));
        %     res_eq(k,:) = double(vpa(subs([Eq1;Eq2;Eq3;Eq4],[S,D,N,Su],sol_temp(ind_m,:))));
        x_eq = sol_temp(ind_m,:);
  
    diff = S_eq-x_eq(1);
    diff_perc = abs(diff)./S_eq;

    if max(diff_perc)>VAR.threshold
        pen_int = 1000*sum(diff(diff_perc>VAR.threshold)); %1000*
    else
        pen_int = 0;
%         save results_pso_generic x_eq lam
    end

    if VAR.k<VAR.N_shell && lam_S<VAR.lambda_prev
        pen_lambda = 10000*(VAR.lambda_prev-lam_S);
    else
        pen_lambda = 0;
    end


    J = VAR.beta1/x_eq(1)+pen_int+pen_lambda; 
    % J = 1e9/sum(lam); % it works but sometimes it increases the number of debris at equilibrium

catch

    J = 1e14;

end


if nargout==3
    varargout{1} = x_eq;
    varargout{2} = diff_perc;
end


end

