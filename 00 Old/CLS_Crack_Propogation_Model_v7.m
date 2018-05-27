%% Matlab reset

clear all
close all
clc;

%% User input

% Support boundary conditions
BC = 'RollerRoller';
%BC = 'ClampedClamped';

% Discretize option
DiscretizeMethod = 'LeftBoundary';
%DiscretizeMethod = 'Central';
%DiscretizeMethod = 'RightBoundary';

% Maximum number of load cycles
N = 10^6; % [-]

% Adherent
E   = 70e9; % [Pa]
t   = 0.0032; % [m]
v   = 0.34; % [-]
% Strap (upper, long adherent) length
L_1 = 0.25; % [m]
% lap (lower, short adherent) length
L_2 = 0.2; % [m]

% Adhesive
t_a = 0.0003; % [m]
E_a = 3.1e9; % [Pa]
G_a = 1.1e9; % [Pa]
v_a = 0.4; % [-]

% Load
S_max = 250e3/t; % [N/m^2]
R_load = 0.1; % [-]

% Initial crack length
b = 0; % [m]

% Number of adherent elements
q = 5000;

%% Load parameters

% Applied force
P_max   = S_max*t;          % [N/m]
P_min   = P_max*R_load;     % [N/m]
P       = [P_min P_max];    % [N/m]

% Applied stress
S_min   = S_max*R_load;     % [N/m^2]
S       = [S_min S_max];    % [N/m^2]

%% Geometric parameters

% Adherent stiffnesses
D_11 = E*t^3/12;
A_11 = E*t;
D_00 = 2*E*(t^3/12+t*(t/2+t_a/2)^2) + E_a*t_a^3/12;
A_00 = E*t*2 + E_a*t_a;

% Initital geometry
l_A0 = L_1-(L_2-b);
l_B0 = L_2-b;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t + t_a)./(2*(l_A0+l_B0));

% Element length
dla = l_A0/q;
dlb = l_B0/q;

%% Discretize adherents
switch DiscretizeMethod
    case 'LeftBoundary'
        X1 = 0:dla:l_A0-dla;
        X0 = 0:dlb:l_B0-dlb;
    case 'Central'
        X1 = dla/2:dla:l_A0-dla/2;
        X0 = dlb/2:dlb:l_B0-dlb/2;
    case 'RightBoundary'
        X1 = dla:dla:l_A0;
        X0 = dlb:dlb:l_B0;
end

%% Calculate (1) overlap edge loads and (2) adhesive stresses
cnt = 0;
for i = 1:q
    
    % Update x0 vector
    x0 = X0(i:end) - (i-1)*dlb;
    
    % Free adherent (l_A) and overlap (l_B) length
    l_A = l_A0+i*dlb;
    l_B = l_B0-i*dlb;
    
    % Overlap edge loads (minumum and maximum)
    [M_kmin, Q_kmin] = Overlap_Edge_Loads(x0, P(1), D_11, D_00, l_A, l_B, t, t_a, BC);
    [M_kmax, Q_kmax] = Overlap_Edge_Loads(x0, P(2), D_11, D_00, l_A, l_B, t, t_a, BC);
    
%     % Adhesive Stresses
     x00 = x0-l_B;
     F = P*cos(alpha);
     [Shear_amin, Peel_amin] = Adhesive_Stresses_LuoTong(x00, F(1), M_kmin, Q_kmin, D_11, A_11, D_00, l_B, E, t, E_a, G_a, t_a, BC);
     [Shear_amax, Peel_amax] = Adhesive_Stresses_LuoTong(x00, F(2), M_kmax, Q_kmax, D_11, A_11, D_00, l_B, E, t, E_a, G_a, t_a, BC);
%     
%     %% Crack growth
%     
%     % SERR
%     G_Imin  = t_a/(2*E_a)*max(Peel_amin,[],2).^2;
%     G_IImin = t_a/(2*G_a)*max(Shear_amin,[],2).^2;
%     G_Tmin  = G_Imin+G_IImin;
%     
%     G_Imax  = t_a/(2*E_a)*max(Peel_amax,[],2).^2;
%     G_IImax = t_a/(2*G_a)*max(Shear_amax,[],2).^2;
%     G_Tmax  = G_Imax+G_IImax;
%     
%     % Mode ratio
%     MR_min  = G_IImin./G_Tmin;
%     MR_max  = G_IImax./G_Tmax;
%     
%     % Disbond growth rate based an D. Burger (2005)
%     % Equivalent mode I
%     G1_eqmin = sqrt(G_Imin)/2+sqrt(G_Imin/2+G_IImin);
%     G1_eqmax = sqrt(G_Imax)/2+sqrt(G_Imax/2+G_IImax);
%     
%     dG1_eq = (G1_eqmax-G1_eqmin).^2;
%     
%     % Paris Law fit coefficients (D. Burger data on FM94)
%     % TO DO: (1) include threshold and (2) G_c
%     c0 = 5.27*10^(-17);
%     m0 = 3.78;
%     c100 = 10^(-17.6);
%     
%     % Note: MR different at P_min and _P_max due to geometric non linearity
%     dbdN = c100.^MR_max.*c0.^(1-MR_max).*dG1_eq.^m0;
end

