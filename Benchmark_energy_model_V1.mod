//Code corresponding Benchmark_energy_model (25 July 2025)

//----------------------------------------------------------------
// Defining Variables and Parameters
//----------------------------------------------------------------

//Variables
var
    C           ${C}$                       (long_name = 'Consumption')
    I           ${I}$                       (long_name = 'Investment')
    r           ${r}$                       (long_name = 'real interest rate')
    K           ${K}$                       (long_name = 'Capital')
    Pi          ${\pi}$                     (long_name = 'Inflation')
    W           ${W}$                       (long_name = 'Wage')
    L           ${L}$                       (long_name = 'Labour')
    Psi         ${\Psi}$                    (long_name = 'Intermediate good firm profit')
    T           ${T}$                       (long_name = 'Tax')
    Omega_K     ${\Omega^K}$                (long_name = 'Investment adjustment cost')
    lambda      ${\lambda}$                 (long_name = 'Lagrange multiplier')
    lambda_K    ${\lambda^K}$               (long_name = 'Lagrange multiplier for capital equation of motion')
    Chi         ${\chi}$                    (long_name = 'Energy production')    
    K_Chi       ${K^\chi}$                  (long_name = 'Energy capital')
    M           ${M}$                       (long_name = 'Stock of pollution')
    Z           ${Z}$                       (long_name = 'Current emissions')
    U           ${U}$                       (long_name = 'Abatement effort')
    Xi          ${\Xi}$                     (long_name = 'Enforced Abatement function')
    I_Chi       ${I^\chi}$                  (long_name = 'Energy investment')
    r_Chi       ${r^\chi}$                  (long_name = 'Rate of return to energy capital')
    Pi_Chi      ${\pi^\chi}$                (long_name = 'Energy sector inflation')
    Tau_Z       ${\tau_Z}$                  (long_name = 'Lagrange multiplier')
    Psi_Chi     ${\Psi^\chi}$               (long_name = 'Energy firm profit')
    MC_Chi      ${MC^\chi}$                 (long_name = 'Energy firm marginal cost')    
    Omega_Chi   ${\Omega^\chi}$             (long_name = 'Energy firm quadratic adjustment cost')    
    Y           ${Y}$                       (long_name = 'Output')
    A           ${A}$                       (long_name = 'Technology shock')
    P_Chi       ${P^\chi}$                  (long_name = 'Lagrange multiplier')
    MC          ${MC}$                      (long_name = 'Intermediate good firm marginal cost')
    Omega       ${\Omega}$                  (long_name = 'Intermediate good firm quadratic adjustment cost ')
    G           ${G}$                       (long_name = 'Government spending')
    Omega_Kchi  ${\Omega^{K^\chi}}$         (long_name = 'Energy investment adjustment cost')
    C_A         ${C_A}$                     (long_name = 'Cost of abatement')
    Eps_IChi    ${\epsilon_{I^\chi}}$
    C_Y
    I_K
    I_Chi_K_Chi
    K_Y
    K_Chi_Chi
    P_ChiChi_Y
;   
//Exogenous variables

varexo 
    Eps_A
    Nu_IChi
    Eps_Z
    ;

//Parameters
    
parameters
    delta
    h
    beta
    omega_K
    gamma
    gamma0
    gamma1
    gamma2
    eta1
    eta2
    xi
    varphi
    phi1
    phi2
    delta_M
    omega_Chi
    pi_Chi
    epsilon 
    alpha_L
    alpha_K
    omega
    pi
    omega_Kchi
    g
    r_ss
    theta_pi
    rho_r
    rho_A
    A_ss
    rho_IChi
    rho_tau
    Tau_Zss
    omega_NKPC
    alpha_Chi
    delta_Chi
    Eps_IChiss
    I_Chiss
    theta_Z
    Z_star
;

//----------------------------------------------------------------
// Calibration
//----------------------------------------------------------------
  
    delta                   = 0.025;         
    h                       = 0.7;              
    beta                    = 0.99;          
    omega_K                 = 4; 
    gamma                   = 1.197;
    gamma0                  = 0.00295;
    gamma1                  = -0.0000066722;
    gamma2                  = 0.000000014647;
    eta1                    = 0.2;          
    eta2                    = 0.1;          
    xi                      = 0.2;            
    varphi                  = 0.45;
    phi1                    = 0.13; 
    phi2                    = 2.6; 
    delta_M                 = 0.025;%0.0035;
    omega_Chi               = 8; 
    pi_Chi                  = 1;      
    epsilon                 = 2;
    alpha_L                 = 0.4;
    alpha_K                 = 0.27;      
    omega                   = 100; 
    pi                      = 1;         
    omega_Kchi              = 8; 
    g                       = 0.13; 
    r_ss                    = 0.035101;
    theta_pi                = 1.5;
    rho_A                   = 0.9;        
    A_ss                    = 1;           
    rho_IChi                = 0.8;     
    rho_tau                 = 0.8;
    Tau_Zss                 = 0.04;
    omega_NKPC              = 6;
    alpha_Chi               = 0.33;
    delta_Chi               = 0.02;
    Eps_IChiss              = 1;
    I_Chiss                 = 0.3;
    theta_Z                 = -1.5;
    Z_star                  = 0.3;
    rho_r                   = 0.6;
    

//----------------------------------------------------------------
// Model
//----------------------------------------------------------------

%%Equations in Model Block Labeled to match the Appendix in the PDF

model;
    
%% Households

    C + I = (r(-1) * K(-1))/Pi + W * L + Psi - T; %(1)

    K = (1 - delta) * K(-1) + (1 - Omega_K) * I; %(2)

    lambda = (C - h * C(-1))^-1 - beta * h * (C(+1) - h * C)^-1; %(3)
	
    lambda_K = beta * (lambda(+1) * r/Pi(+1) + lambda_K(+1) * (1 - delta)); %(4)

    lambda = lambda_K * (1 - Omega_K - (omega_K/I(-1)) * ((I/I(-1)) - 1) * I) + beta * lambda_K(+1) * omega_K * ((I(+1)/I)^3 - (I(+1)/I)^2); %(5)

    lambda * W = gamma * (1 - L)^(-1); %(6)   

    Omega_K = (omega_K/2) * (I/I(-1) - 1)^2; %(7)

	% Eqs (8) & (9) substituted in

%% Energy producer

	Chi = A* (1 - (gamma0 + gamma1*M + gamma2*M^2)) * K_Chi(-1)^alpha_Chi; %(10)

    % Eq (11) substituted in

    Z = (1 - (U + Xi)) * varphi * Chi; %(12)

    Xi = (eta1 * log(1 + xi*I_Chi))/(1 + eta2 * log(1 + xi*I_Chi)); %(13)
	
    C_A = phi1 * (U + Xi)^phi2 * Chi; %(14)

    M = (1 - delta_M) * M(-1) + Z; %(15)

    r_Chi/Pi_Chi(+1) = (MC_Chi(+1)) * alpha_Chi * (Chi(+1))/K_Chi; %(16)
	
	varphi * Tau_Z = phi1 * phi2 * (U + Xi)^(phi2 - 1); %(17)

    Psi_Chi = P_Chi * Chi - r_Chi(-1)*K_Chi(-1)/Pi_Chi - C_A - Tau_Z*Z; %(18)

	MC_Chi = (((1 - (gamma0 + gamma1*M + gamma2*M^2)) * A)^(-1)) * ((r_Chi(-1)/(alpha_Chi*Pi_Chi))^alpha_Chi) + phi1 * (U + Xi)^phi2 + Tau_Z * (1 - (U + Xi)) * varphi; %(19)

	Omega_Chi = (omega_Chi/2) * (Pi_Chi/pi_Chi - 1)^2; %(20)

    omega_Chi * (Pi_Chi/pi_Chi - 1) * (Pi_Chi/pi_Chi) = beta * lambda(+1)/lambda * omega_Chi * (Pi_Chi(+1)/pi_Chi - 1) * (Pi_Chi(+1)/pi_Chi) * (Chi(+1)/Chi) + epsilon*MC_Chi + (1 - epsilon); %(21) 
	
%% Intermediate good firms

    Y = (1 - (gamma0 + gamma1*M + gamma2*M^2)) * A * (L^alpha_L) * (K(-1)^alpha_K) * (Chi^(1 - alpha_L - alpha_K)); %(22)

    r/Pi(+1) = MC(+1) * alpha_K * (Y(+1)/K); %(23)
	
    W = MC * alpha_L * (Y/L); %(24)	
	
    P_Chi = MC * (1 - alpha_L - alpha_K) * (Y/Chi); %(25)
	
    Psi = Y - W*L - (r(-1)*K(-1))/Pi - P_Chi * Chi; %(26)
	
    Omega = (omega/2) * (Pi/pi - 1)^2; %(27)

    omega * (Pi/pi - 1) * (Pi/pi) = beta * lambda(+1)/lambda * omega * (Pi(+1)/pi - 1) * (Pi(+1)/pi) * (Y(+1)/Y) + omega_NKPC*MC + (1 - omega_NKPC); %(28)
   	
%% Government 
    
    G + I_Chi = T + (r_Chi(-1) * K_Chi(-1))/Pi_Chi + Tau_Z*Z + Psi_Chi; %(29)

    G = g * Y; %(30)

    K_Chi = (1 - delta_Chi) * K_Chi(-1) + (1 - Omega_Kchi)*(1 - xi)*I_Chi; %(31)

    Omega_Kchi = (omega_Kchi/2) * (I_Chi/I_Chi(-1) - 1)^2; %(32)

    I_Chi/I_Chiss  = (Z/Z_star)^(theta_Z) * Eps_IChi; %(33)

    (1+r)/(1+r_ss) = (Pi/pi)^theta_pi; %(34)

	
%% Model shocks

    log(A) = (1 - rho_A) * log(A_ss) + rho_A * log(A(-1)) + Eps_A; %(35)

    log(Eps_IChi) = (1-rho_IChi) * log(Eps_IChiss) + rho_IChi * log(Eps_IChi(-1)) + Nu_IChi; %(36)

    log(Tau_Z) = (1 - rho_tau) * log(Tau_Zss) + rho_tau * log(Tau_Z(-1)) + Eps_Z; %(37)

%% SSs

    C_Y = C/Y;

    I_K = I/K;

    I_Chi_K_Chi = I_Chi/K_Chi;

    K_Y = K/Y;

    K_Chi_Chi = K_Chi/Chi;

    P_ChiChi_Y = (P_Chi*Chi)/Y;

end;

model_diagnostics;
resid;
steady;
check;


//----------------------------------------------------------------
// Shocks
//----------------------------------------------------------------

 shocks;
     var Eps_A;     stderr 0.01; // Technology shock
     var Eps_Z;     stderr 0.01; // Tax on energy firm emissions shock 
     var Nu_IChi;   stderr 0.01; // Government investment in energy firms shock 
 end;

//----------------------------------------------------------------
// Solving the Model
//----------------------------------------------------------------

%write_latex_original_model;
stoch_simul(IRF = 40); 


% %----------------------------------------------------------------------%
% % 0. GLOBAL FIGURE SETTINGS – one-time, no functions, no size bugs
% %----------------------------------------------------------------------%
% irf_horizon = options_.irf;                         % length of IRFs
% 
% set(0,'DefaultFigureUnits','centimeters');          % <- all sizes are cm
% set(0,'DefaultFigurePosition',[2 2 16 12]);         %   (x y w h)
% set(0,'DefaultFigureColor','w');
% set(0,'DefaultFigurePaperUnits','centimeters');
% set(0,'DefaultFigurePaperPositionMode','auto');     % on-screen = on-paper
% 
% set(0,'defaultAxesTickLabelInterpreter','latex');
% set(0,'defaultTextInterpreter','latex');
% set(0,'defaultLegendInterpreter','latex');
% 
% %----------------------------------------------------------------------%
% % 1. LaTeX variable labels (containers.Map so we can call by name)
% %----------------------------------------------------------------------%
% VARIABLE_LABELS = containers.Map( ...
%  {'Y','C','I','L','Pi','r','K','W','Psi','T', ...
%   'Omega_K','lambda','lambda_K','Chi','A_Chi','K_Chi', ...
%   'M','Z','U','Xi','I_Chi','r_Chi','Pi_Chi','Tau_Z', ...
%   'Psi_Chi','MC_Chi','Omega_Chi','P_Chi','MC','Omega', ...
%   'G','Omega_Kchi','C_A'}, ...
%  {'Y_t','C_t','I_t','L_t','\pi_t','r_t','K_t','W_t', ...
%   '\Psi_t','T_t','\Omega_t^K','\lambda_t','\lambda_t^K', ...
%   '\chi_t','A_t^\chi','K_t^\chi','M_t','Z_t','U_t','\Xi_t', ...
%   'I_t^\chi','r_t^\chi','\pi_t^\chi','\tau_{Z,t}','\Psi_t^\chi', ...
%   'MC_t^\chi','\Omega_t^\chi','P_t^\chi','MC_t','\Omega_t', ...
%   'G_t','\Omega_t^{K^\chi}','C_{A,t}'});
% 
% %----------------------------------------------------------------------%
% % 2.  SPECIFY (i) variable lists, (ii) shock codes, (iii) file names
% %----------------------------------------------------------------------%
% specs = { ...
%   struct('vars',{{'Y','C','I','L','r','r_Chi','Pi','Pi_Chi', ...
%                   'Chi','C_A','MC','MC_Chi','M','Z','Xi','U'}}, ...
%          'shock','Eps_A',   ...
%          'title','IRFs to Technology Shock ($\epsilon_A$)', ...
%          'file','IRF_TechShock'), ...
%   struct('vars',{{'Y','C','I','L','r','r_Chi','Pi','Pi_Chi', ...
%                   'Chi','C_A','MC','MC_Chi','M','Z','Xi','U'}}, ...
%          'shock','Nu_IChi', ...
%          'title','IRFs to Energy Investment Shock ($\nu_{I^\chi}$)', ...
%          'file','IRF_EnergyInvestmentShock'), ...
%   struct('vars',{{'Y','C','I','L','r','r_Chi','Pi','Pi_Chi', ...
%                   'Chi','C_A','MC','MC_Chi','M','Z','Xi','U'}}, ...
%          'shock','Eps_Z',   ...
%          'title','IRFs to Emission Tax Shock ($\Eps_{\tau_Z}$)', ...
%          'file','IRF_EmissionTaxShock')};
% 
% %----------------------------------------------------------------------%
% % 3.  LOOP THROUGH THE FOUR FIGURES
% %----------------------------------------------------------------------%
% for s = 1:numel(specs)
%     V   = specs{s}.vars;
%     shk = specs{s}.shock;
%     ttl = specs{s}.title;
%     out = specs{s}.file;
% 
%     figure('Name',ttl);                         % size is 16 × 12 cm
%     for k = 1:numel(V)
%         subplot(4,4,k);
%         v   = V{k};
%         fld = [v '_' shk];
%         if ~isfield(oo_.irfs,fld)
%             title(['$' VARIABLE_LABELS(v) '$ (no data)'],'Interpreter','latex');
%             axis off; continue;
%         end
%         plot(1:irf_horizon,oo_.irfs.(fld),'-k','LineWidth',1.5);
%         title(['$',VARIABLE_LABELS(v),'$'],'Interpreter','latex');
%         xlabel('Periods','Interpreter','latex'); grid on;
%     end
% 
%     % ------- export as VECTOR EPS with correct bounding box ------------
%     print('-depsc2','-loose',[out '.eps']);     % EPS only, no -r needed
%     %close(gcf);                                % tidy up (optional)
% end