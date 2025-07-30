% Benchmark_energy_model_V1_steadystate.m  - must sit in the same directory as Benchmark_energy_model_V1.mod
function [ys,params,check] = Benchmark_energy_model_V1_steadystate(ys,~,M_,~)

%% ---------- 0. boiler-plate ----------
check  = 0;
params = M_.params;                     % local copy of parameters

%Giving the parameters from Dynare
delta       = params(strcmp('delta',       M_.param_names));
h           = params(strcmp('h',           M_.param_names));
beta        = params(strcmp('beta',        M_.param_names));
omega       = params(strcmp('omega',       M_.param_names));
omega_K     = params(strcmp('omega_K',     M_.param_names));
omega_Chi   = params(strcmp('omega_Chi',   M_.param_names));
omega_Kchi  = params(strcmp('omega_Kchi',  M_.param_names));
omega_NKPC  = params(strcmp('omega_NKPC',  M_.param_names));
gamma       = params(strcmp('gamma',       M_.param_names));
gamma0      = params(strcmp('gamma0',      M_.param_names));
gamma1      = params(strcmp('gamma1',      M_.param_names));
gamma2      = params(strcmp('gamma2',      M_.param_names));
eta1        = params(strcmp('eta1',        M_.param_names));
eta2        = params(strcmp('eta2',        M_.param_names));
xi          = params(strcmp('xi',          M_.param_names));
varphi      = params(strcmp('varphi',      M_.param_names));
phi1        = params(strcmp('phi1',        M_.param_names));
phi2        = params(strcmp('phi2',        M_.param_names));
delta_M     = params(strcmp('delta_M',     M_.param_names));
delta_Chi   = params(strcmp('delta_Chi',   M_.param_names));
pi          = params(strcmp('pi',          M_.param_names));
pi_Chi      = params(strcmp('pi_Chi',      M_.param_names));
epsilon     = params(strcmp('epsilon',     M_.param_names));
alpha_L     = params(strcmp('alpha_L',     M_.param_names));
alpha_K     = params(strcmp('alpha_K',     M_.param_names));
alpha_Chi   = params(strcmp('alpha_Chi',   M_.param_names));
g           = params(strcmp('g',           M_.param_names));
r_ss        = params(strcmp('r_ss',        M_.param_names)); % Does not need to be initialised in the MATLAB code but is here for completeness
rho_pi      = params(strcmp('rho_pi',      M_.param_names)); % Does not need to be initialised in the MATLAB code but is here for completeness
rho_A       = params(strcmp('rho_A',       M_.param_names)); % Does not need to be initialised in the MATLAB code but is here for completeness
rho_IChi    = params(strcmp('rho_IChi',    M_.param_names)); % Does not need to be initialised in the MATLAB code but is here for completeness
rho_tau     = params(strcmp('rho_tau',     M_.param_names)); % Does not need to be initialised in the MATLAB code but is here for completeness
A_ss        = params(strcmp('A_ss',        M_.param_names));
Tau_Zss     = params(strcmp('Tau_Zss',     M_.param_names));
Eps_IChiss  = params(strcmp('Eps_IChiss',  M_.param_names));
I_Chiss     = params(strcmp('I_Chiss',     M_.param_names)); 
theta_Z     = params(strcmp('theta_Z',     M_.param_names)); 
Z_star      = params(strcmp('Z_star',      M_.param_names));
rho_r       = params(strcmp('rho_r',       M_.param_names)); % Does not need to be initialised in the MATLAB code but is here for completeness

MC_Chiss = (epsilon-1)/epsilon;

%% ---------- 1. numerical core via fsolve ----------

% x = [lambda_K , I , L, K_Chi; I_Chi; M; U; Xi]     % 8 unknowns
x0 = [2 ; 0.15 ; 0.33; 5; 0.15; 50; 0.3; 0.001];     % rough guess

opt = optimoptions('fsolve','Display','off','TolFun',1e-14,'TolX',1e-14);

% note: core_residual is a **nested function** defined just below
[x_ss,~,exitflag] = fsolve(@core_residual,x0,opt);

if exitflag<=0
    check = 1;                          % tell Dynare steady state failed
    return
end
lambda_K = x_ss(1);   I = x_ss(2);   L = x_ss(3); K_Chi = x_ss(4); I_Chi = x_ss(5); M = x_ss(6); U = x_ss(7); Xi = x_ss(8);

%% ---------- 2. rebuild the full steady state ----------

%Giving the steady-state relationships from the model  

lambda_K = x_ss(1);
I        = x_ss(2);
L        = x_ss(3);
K_Chi    = x_ss(4);
I_Chi    = x_ss(5); 
M        = x_ss(6);
U        = x_ss(7); 
Xi       = x_ss(8);
Z = Z_star;
Eps_IChi = Eps_IChiss;
lambda = lambda_K;
C = (1 - beta*h)/(lambda*(1-h));
K = I/delta; 
A = A_ss;
Tau_Z = Tau_Zss;
Pi_Chi = pi_Chi;
Pi = pi;
MC = (omega_NKPC - 1)/omega_NKPC;
MC_Chi = MC_Chiss;
Omega_Chi = (omega_Chi/2) * ((Pi_Chi/pi_Chi - 1)^2);
Omega = (omega/2) * (Pi/pi - 1)^2;
Chi = A * (1 - (gamma0 + gamma1*M + gamma2*M^2)) * K_Chi^alpha_Chi;
Z = (1 - (U + Xi)) * varphi * Chi;
C_A = phi1 * (U + Xi)^phi2 * Chi;
Y = (1 - (gamma0 + gamma1*M + gamma2*M^2)) * A * (L^alpha_L) * (K^alpha_K) * (Chi^(1 - alpha_L - alpha_K));
P_Chi = MC * (1 - alpha_L - alpha_K) * (Y/Chi);
r_Chi = Pi_Chi * (MC_Chi) * alpha_Chi * (Chi)/K_Chi;
W = MC * alpha_L * (Y/L);
r = Pi * (1/beta - (1 - delta));
Psi = Y - W*L - (r*K)/Pi - P_Chi * Chi;
Psi_Chi = (P_Chi * Chi - r_Chi*K_Chi/Pi_Chi - C_A - Tau_Z*Z);
G = g * Y;
T = G + I_Chi - (r_Chi * K_Chi)/Pi_Chi - Tau_Z * Z - Psi_Chi;
Omega_Kchi = (omega_Kchi/2) * (I_Chi/I_Chi - 1)^2;
Omega_K = (omega_K/2) * (I/I - 1)^2;
C_Y = C/Y;
I_K = I/K;
I_Chi_K_Chi = I_Chi/K_Chi;
K_Y = K/Y;
K_Chi_Chi = K_Chi/Chi;
P_ChiChi_Y = (P_Chi*Chi)/Y;

I_Chiss = I_Chi / ((Z/Z_star)^theta_Z * Eps_IChi);     % 1. overwrite local copy
params(strcmp('I_Chiss',M_.param_names)) = I_Chiss;    % 2. store in param list
ys(strcmp('I_Chiss',M_.endo_names)) = I_Chiss;         % 3. store in state vector

% I_Chiss = (delta_Chi * K_Chi)/(1-xi);                  % 1. overwrite local copy
% params(strcmp('I_Chiss',M_.param_names)) = I_Chiss;    % 2. store in param list
% ys(strcmp('I_Chiss',M_.endo_names)) = I_Chiss;         % 3. store in state vector

%% ---------- 3. load results into ys ----------
ys(strcmp('lambda_K',     M_.endo_names)) = lambda_K;
ys(strcmp('I',            M_.endo_names)) = I;
ys(strcmp('L',            M_.endo_names)) = L;
ys(strcmp('M',            M_.endo_names)) = M;
ys(strcmp('lambda',       M_.endo_names)) = lambda;
ys(strcmp('C',            M_.endo_names)) = C;
ys(strcmp('K',            M_.endo_names)) = K;
ys(strcmp('A',            M_.endo_names)) = A;
ys(strcmp('Tau_Z',        M_.endo_names)) = Tau_Z;
ys(strcmp('Pi_Chi',       M_.endo_names)) = Pi_Chi;
ys(strcmp('Pi',           M_.endo_names)) = Pi;
ys(strcmp('MC',           M_.endo_names)) = MC;
ys(strcmp('MC_Chi',       M_.endo_names)) = MC_Chi;
ys(strcmp('Omega_Chi',    M_.endo_names)) = Omega_Chi;
ys(strcmp('Omega',        M_.endo_names)) = Omega;
ys(strcmp('I_Chi',        M_.endo_names)) = I_Chi;
ys(strcmp('K_Chi',        M_.endo_names)) = K_Chi;
ys(strcmp('Chi',          M_.endo_names)) = Chi;
ys(strcmp('Xi',           M_.endo_names)) = Xi;
ys(strcmp('U',            M_.endo_names)) = U;
ys(strcmp('Z',            M_.endo_names)) = Z;
ys(strcmp('r_Chi',        M_.endo_names)) = r_Chi;
ys(strcmp('C_A',          M_.endo_names)) = C_A;
ys(strcmp('Psi_Chi',      M_.endo_names)) = Psi_Chi;
ys(strcmp('Y',            M_.endo_names)) = Y;
ys(strcmp('W',            M_.endo_names)) = W;
ys(strcmp('r',            M_.endo_names)) = r;
ys(strcmp('P_Chi',        M_.endo_names)) = P_Chi;
ys(strcmp('Psi',          M_.endo_names)) = Psi;
ys(strcmp('G',            M_.endo_names)) = G;
ys(strcmp('T',            M_.endo_names)) = T;
ys(strcmp('Omega_Kchi',   M_.endo_names)) = Omega_Kchi;
ys(strcmp('Omega_K',      M_.endo_names)) = Omega_K;
ys(strcmp('C_Y',          M_.endo_names)) = C_Y;
ys(strcmp('I_K',          M_.endo_names)) = I_K;
ys(strcmp('I_Chi_K_Chi',  M_.endo_names)) = I_Chi_K_Chi;
ys(strcmp('K_Y',          M_.endo_names)) = K_Y;
ys(strcmp('K_Chi_Chi',    M_.endo_names)) = K_Chi_Chi;
ys(strcmp('P_ChiChi_Y',   M_.endo_names)) = P_ChiChi_Y;
ys(strcmp('Eps_IChi',     M_.endo_names)) = Eps_IChi;

%% ========== NESTED RESIDUAL FUNCTION =================================
    function F = core_residual(x)

    lambda_K = x(1);  I = x(2);  L = x(3); K_Chi = x(4); I_Chi = x(5); M = x(6); U = x(7); Xi = x(8);

        % ---- reuse algebra that only depends on x & params ----
        Eps_IChi = Eps_IChiss;
        lambda = lambda_K;
        C = (1 - beta*h)/(lambda*(1-h));
        K = I/delta; 
        A = A_ss;
        Tau_Z = Tau_Zss;
        Pi_Chi = pi_Chi;
        Pi = pi;
        Z = Z_star;
        MC = (omega_NKPC - 1)/omega_NKPC;
        MC_Chi = MC_Chiss;
        Omega_Chi = (omega_Chi/2) * ((Pi_Chi/pi_Chi - 1)^2);
        Omega = (omega/2) * (Pi/pi - 1)^2;
        Chi = A * (1 - (gamma0 + gamma1*M + gamma2*M^2)) * K_Chi^alpha_Chi;
        Z = (1 - (U + Xi)) * varphi * Chi;
        C_A = phi1 * (U + Xi)^phi2 * Chi;
        Y = (1 - (gamma0 + gamma1*M + gamma2*M^2)) * A * (L^alpha_L) * (K^alpha_K) * (Chi^(1 - alpha_L - alpha_K));
        P_Chi = MC * (1 - alpha_L - alpha_K) * (Y/Chi);
        r_Chi = Pi_Chi * (MC_Chi) * alpha_Chi * (Chi)/K_Chi;
        W = MC * alpha_L * (Y/L);
        r = Pi * (1/beta - (1 - delta));
        Psi = Y - W*L - (r*K)/Pi - P_Chi * Chi;
        Psi_Chi = (P_Chi * Chi - r_Chi*K_Chi/Pi_Chi - C_A - Tau_Z*Z);
        G = g * Y;
        T = G + I_Chi - (r_Chi * K_Chi)/Pi_Chi - Tau_Z * Z - Psi_Chi;
        Omega_Kchi = (omega_Kchi/2) * (I_Chi/I_Chi - 1)^2;
        Omega_K = (omega_K/2) * (I/I - 1)^2;
        C_Y = C/Y;
        I_K = I/K;
        I_Chi_K_Chi = I_Chi/K_Chi;
        K_Y = K/Y;
        K_Chi_Chi = K_Chi/Chi;
        P_ChiChi_Y = (P_Chi*Chi)/Y;

        % ---- 4 residual equations  ------
        F          = zeros(8,1);
        F(1)       = I - ((r * K)/Pi + W * L + Psi - T - C); % solving for lambda_K
        F(2)       = r/Pi - MC * alpha_K * (Y/K);   % For solving I
        F(3)       = L - (1 - gamma/(lambda * W)); % For solving L
        F(4)       = MC_Chi - ((((1 - (gamma0 + gamma1*M + gamma2*M^2)) * A)^(-1)) * ((r_Chi/(alpha_Chi*Pi_Chi))^alpha_Chi) + phi1 * (U + Xi)^phi2 + Tau_Z * (1 - (U + Xi)) * varphi);
        F(5)       = I_Chi - (delta_Chi * K_Chi)/(1-xi);
        F(6)       = M - Z/delta_M;
        F(7)       = Xi - (eta1 * log(1 + xi*I_Chi))/(1 + eta2 * log(1 + xi*I_Chi));
        F(8)       = U - (((varphi * Tau_Z)/(phi1 * phi2))^(1/(phi2-1)) - Xi);
    end
% ========== end of nested function ====================================
end