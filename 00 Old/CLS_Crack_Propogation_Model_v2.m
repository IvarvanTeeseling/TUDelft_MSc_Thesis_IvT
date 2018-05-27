clear all; close all; clc

%% Input

% Load
S_max = 300e6; % [N/m^2]
R_load = 0.1; % [-]

% Total number of cycles
N = 10; % [-]

% Adherent
E   = 70e9; % [Pa]
t   = 0.0032; % [m]
v   = 0.34; % [-]
% Strap (upper, long adherent) length
L_1 = 1; % [m]
% lap (lower, short adherent) length
L_2 = 0.5; % [m]

% Adhesive
t_a = 0.078*0.0032; % [m]
E_a = 0.04*70e9; % [Pa]
G_a = 1.1e9; % [Pa]
v_a = 0.4; % [-]

% Initial crack length
b = 0; % [m]

%% Calculations

% Applied force
P_max = S_max*t; % [N/m]
P_min = P_max*R_load; % [N/m]
P_mean = P_max*(1+R_load)/2; % [N/m]
P_load = [P_min P_max P_mean]';

% Applied stress
S_min = S_max*R_load; % [N/m^2]
S_mean = S_max*(1+R_load)/2; % [N/m^2]
S_load = [S_min S_max S_mean];

% Adherent stiffnesses
% Free adherent
D_11 = E*t^3/12;
A_11 = E*t;
% Overlap region (lumped together)
D_00 = 2*E*(t^3/12+t*(t/2+t_a/2)^2) + E_a*t_a^3/12;
A_00 = E*t*2 + E_a*t_a;

% Calculate displacement and overlap edge loads
G_I     = zeros(1,2);
G_II    = zeros(1,2);
G_T     = zeros(1,2);
G1_eq   = zeros(1,2);
MR      = zeros(1,2);
N = 0:0.01:0.5;
b_hist  = zeros(1,length(N));
G_hist  = zeros(4,length(N));

% Initital geometry
l_A = L_1-(L_2+b);
l_B = L_2-b;

% Initial angle neutral axis w.r.t. adherent neutral axis
alpha = (t + t_a)/(2*(l_A+l_B));

% x-vector spanning the free adherent length
x1 = 0:l_A/100:l_A;
% x-vector spanning the overlap length
x0 = 0:l_B/100:l_B;
% x-vector spanning the entire length
x = [x1 x0+l_A];

for i = 1:2
    
    P = P_load(i);
    
    % >> Roller-Roller Boundary Conditions
    
    % Eigenvalues
    lambda_1 = sqrt(P/D_11);
    lambda_0 = sqrt(P/D_00);
    
    % Integration coeficients
    A_1 = 0;
    %B_1 = -sqrt(D_11) * (cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) * t + cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) * t_a + (2 * alpha * l_A) + 0.2e1 * alpha * l_B - t - t_a) / (sqrt(D_11) * cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) * sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) + sqrt(D_00) * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B)) / 0.2e1;
    B_1 = -(sqrt(D_11)*(cosh(lambda_0*l_B)*t+cosh(lambda_0*l_B)*t_a+2*alpha*(l_A+l_B)-t-t_a))/(2*(sqrt(D_11)*cosh(lambda_0*l_B)*sinh(lambda_1*l_A)+sqrt(D_00)*cosh(lambda_1*l_A)*sinh(lambda_0*l_B)));
    %A_0 = (-2*sqrt(D_00/D_11)*B_1*cosh(lambda_1*l_A)*sinh(lambda_0*l_B)-2*alpha*(l_A+l_B)+t+t_a)/(2*cosh(lambda_0*l_B))
    A_0 = (-0.2e1 * sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) * D_11 ^ (-0.1e1 / 0.2e1) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) + (-0.2e1 * l_A - 0.2e1 * l_B) * alpha + t + t_a) / cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) / 0.2e1;
    B_0 = sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) * D_11 ^ (-0.1e1 / 0.2e1);
    B_0 = sqrt(D_00/D_11)*B_1*cosh(lambda_1*l_A)
    
    % From force equilibrium (M_A = M_B = 0) due to roller BC
    R_A = 0;
    M_A = 0;
    
    % Vertical displacement of the free adherent
    w1 = A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha + R_A/P) .* x1 + M_A/P;
    % Vertical desiplacement of the overlap region
    w0 = A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + alpha .* (l_A + x0) - (t + t_a)/2 + (l_A + x0) * R_A/P + M_A/P;
    
    % Bending moment of the free adherent
    % displacement as positive, I think it should be +)
    M1 = P * (-A_1 * cosh(lambda_1 .* x1) - B_1 * sinh(lambda_1 .* x1));
    % Bending moment of the overlap region
    M0 = P * (-A_0 * cosh(lambda_0 .* x0) - B_0 * sinh(lambda_0 .* x0));
    
    % Shear force in the free adherent
    Q1 = P * (-A_1*lambda_1*sinh(lambda_1.*x1)-B_1*lambda_1*sinh(lambda_1.*x1));
    % Shear force in the overlap region
    Q2 = P * (-A_0*lambda_0*sinh(lambda_0.*x0)-B_0*lambda_0*sinh(lambda_0.*x0));
    
    % Overlap edge loads
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 31
    M_k = M1(end);
    Q_k = Q1(end);
    
    % Adhesive Stresses
    x0 = -l_B:l_B/1000:0;
    F = P*cos(alpha);
    
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 30
    beta_k = sqrt(F/D_11);
    beta_0 = sqrt(F/D_00);
    
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 33
    alpha_k = A_11 * t^2 / (4 * D_11);
    alpha_a = (1 + alpha_k) / 4;
    
    % Bending moment factor
    k = M_k/(P*(t+t_a)/2);
    % By GR (note: t_a = 0 in for this expression)
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 40
    k = 1/(1+2*sqrt(2)*tanh(beta_k*l_B/(2*sqrt(2))*coth(beta_k*l_A)));
    
    % Peel stress integration constants
    % > From Luo and Tong 2004
    beta_s = sqrt(2)/2*(24*E_a/(E*t^3*t_a))^(1/4);
    beta_t = sqrt(8*G_a/(A_11*t_a));
    beta_a = sqrt(alpha_a*beta_t^2);
    
    % Peel stress integration constants
    % > From Luo and Tong 2004
    B_s1 = 6/(beta_s^3*E*t^3)*(M_k*beta_s*(sinh(beta_s*l_B)*cos(beta_s*l_B)+cosh(beta_s*l_B)*sin(beta_s*l_B))+Q_k*sinh(beta_s*l_B)*sin(beta_s*l_B))/(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));
    B_s4 = 6/(beta_s^3*E*t^3)*(M_k*beta_s*(sinh(beta_s*l_B)*cos(beta_s*l_B)-cosh(beta_s*l_B)*sin(beta_s*l_B))+Q_k*cosh(beta_s*l_B)*cos(beta_s*l_B))/(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));
    
    % Adhesive shear stress (Luo and Tong 2007)
    Shear_a = beta_t*(F*t+6*M_k)*cosh(beta_t*x0)/(8*t*sinh(beta_t*l_B)) + 3*(F*t-2*M_k)/(8*t*l_B);
    % Adheisve peel stress (Luo and Tong 2004)
    Peel_a = 2*E_a/t_a*(B_s1*sinh(beta_s*x0).*sin(beta_s*x0)+B_s4*cosh(beta_s*x0).*cos(beta_s*x0));
    
    G_1 = t_a/(2*E_a)*Peel_a(1)^2;
    G_2 = t_a/(2*G_a)*Shear_a(1)^2;
    G_t = G_1+G_2;
    mr = G_2/G_t;
    
    G_I(i) = G_1;
    G_II(i) = G_2;
    G_T(i) = G_t;
    MR(i) = mr;
    G1_eq(i) = sqrt(G_1)/2+sqrt(G_1/2+G_2);
end

dG1_eq = (G1_eq(2)-G1_eq(1))^2;

c0 = 5.27*10^(-17);
m0 = 3.78;
c100 = 10^(-17.6);
MR = MR(2);

dbdN = c100^MR*c0^(1-MR)*dG1_eq^m0;

figure(2)
hold on
plot(x0/l_B,Shear_a/10^6,'b')
hold off

figure(3)
hold on
plot(x0/l_B,Peel_a/10^6,'b')
hold off

% k1 = (1+t_a/t)*k;
% Shear_a_max = 1/8*((3*k1+1)*beta_t*coth(beta_t*l_B)+3*(1-k1))*F;
% Peel_a_max = (1+(beta_k/beta_s)*coth(beta_k*l_A))*beta_s^2*M_k;

% >> Clamped-Clamped Boundary Conditions