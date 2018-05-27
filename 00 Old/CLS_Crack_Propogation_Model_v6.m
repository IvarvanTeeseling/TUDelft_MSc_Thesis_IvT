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

% Adherent
E   = 70e9; % [Pa]
t   = 0.0016; % [m]
v   = 0.34; % [-]
% Strap (upper, long adherent) length
L_1 = 0.40; % [m]
% lap (lower, short adherent) length
L_2 = 0.10; % [m]

% Adhesive
t_a = 0.000125; % [m]
E_a = 3.1e9; % [Pa]
G_a = 1.1e9; % [Pa]
v_a = 0.4; % [-]

% Load
S_max = 200e3/t; % [N/m^2]
R_load = 0.1; % [-]

% Initial crack length
b = 0; % [m]

% Number of adherent elements
q = 50;

% Max disbond elements
qmax = q * 5/8;
qmax = ceil(qmax);

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
        x1 = 0:dla:l_A0-dla;
        x0 = 0:dlb:l_B0-dlb;
    case 'Central'
        x1 = dla/2:dla:l_A0-dla/2;
        x0 = dlb/2:dlb:l_B0-dlb/2;
    case 'RightBoundary'
        x1 = dla:dla:l_A0;
        x0 = dlb:dlb:l_B0;
end

%% Overlap edge loads

% Create x0 vector for each time an entire element has been cracked
x0 = x0.*ones(qmax,1);
x0 = x0 - dlb*ones(qmax,q).*(0:qmax-1)';
x0(x0<0) = 0;

% Free adherent (l_A) and overlap (l_B) length
l_A = l_A0*ones(qmax,1) + dlb*(0:qmax-1)';
l_B = l_B0*ones(qmax,1) - dlb*(0:qmax-1)';

% Overlap edge loads (minumum and maximum)
[M, Q] = Overlap_Edge_Loads(x0, P, D_11, D_00, l_A, l_B, t, t_a, BC, 3);

% Extract solutions
M_k = M.A;
Q_k = Q.A;

M_0 = M.B;
Q_0 = Q.B;

%% Adhesive Stresses
x00 = x0-l_B;
F = P*cos(alpha);
[Shear_a, Peel_a, Ad_1, Ad_2] = Adhesive_Stresses_LuoTong(x00, F, M_k, Q_k, D_11, A_11, D_00, l_B, E, t, E_a, G_a, t_a, BC);

%% Adherent stresses

% Y-axis pointing down
s1_x0   = Ad_1.N/t;
s1_xb   = Ad_1.M/D_11*t/2;
s1_t    = s1_x0+s1_xb;

s2_x0   = Ad_2.N/t;
s2_xb   = Ad_2.M/D_11*-1*t/2;
s2_t    = s2_x0+s2_xb;

%% Strain Energy Release Rate

% Lai et al. (1995); using the Suo and Hutchinson (1990) approach
% mu = 1;
% sm = 1;
% delta = 1;
% I = sm*((delta-1/mu)^2-(delta-1/mu)+1/3)+delta/mu*(delta-1/mu)+1/(3*mu^3);
% A = 1/mu + sm;
% G_T1 = 1/(2*E)*((P.^2/t+12*M_k.^2./t^3)+(-P.^2/(A*t)-M_0.^2/(I*t^3)));

% Brussat et al (1977); infinitely long specimen 
%G_T2 = P.^2/(4*E*t)

% SERR
G_I  = t_a/(2*E_a)*Shear_a(:,1,:).^2;
G_II = t_a/(2*G_a)*Peel_a(:,1,:).^2;
G_T  = G_I+G_II;

% Mode ratio
% Note: MR different at P_min and _P_max due to geometric non linearity
MR  = G_II./G_T;

%% Disbond growth rate 
% Applying D. Burger (2005) PhD data and model for FM94 adhesive

% Equivalent mode I
G1_eq = sqrt(G_I)/2+sqrt(G_I/2+G_II);

dG1_eq = (G1_eq(:,:,2)-G1_eq(:,:,1)).^2;

% Paris Law fit coefficients (D. Burger data on FM94)
% TO DO: (1) include threshold and (2) G_c
c0 = 5.27*10^(-17);
m0 = 3.78;
c100 = 10^(-17.6);

% Disbond growth rate using the MR at the maximum load
dbdN = c100.^MR(:,:,2).*c0.^(1-MR(:,:,2)).*dG1_eq.^m0;

% Load cycles for the disbond to grow the element distance
n = dlb./dbdN;

%% Fatigue accumulation



%% Plot results 

%plot(x0(1,:),Shear_a(1,:,1))
%plot(cumsum(n),x0(1,1:qmax)+dlb)

% figure(1)
% hold on
% plot(x0(1,:),Ad_1.N(1,:,1),'r')
% plot(x0(1,:),Ad_2.N(1,:,1),'b')
% hold off
% 
% figure(2)
% hold on
% plot(x0(1,:),Ad_1.Q(1,:,1),'r')
% plot(x0(1,:),Ad_2.Q(1,:,1),'b')
% hold off
% 
% figure(3)
% hold on
% plot(x0(1,:),Ad_1.M(1,:,1),'r')
% plot(x0(1,:),Ad_2.M(1,:,1),'b')
% hold off
% 
% figure(4)
% hold on
% plot(x0(1,:),s1_t(1,:,1),'r')
% plot(x0(1,:),s2_t(1,:,1),'b')
% hold off

