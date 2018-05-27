clear all; close all; clc
%% Assumptions
% 1. Symmetric & (quasi) istropic substrates
%   > E, t, v are identical for both adherends
%   > FML laminate must be (quasi) isotropic
% 2. Adhesive thickness t_a is included (contrary to GR)

%% Input

% Applied load along the joint neutral line [N/m]
P   = 20000;

% Strap (upper, long adherent) & lap (lower, short adherent) properties
E   = 72e9; % [Pa]
t   = 0.0032; % [m]
v   = 0.33; % [-]
% Strap length
L_1 = 0.25; % [m]
% Lap length
L_2 = 0.05; % [m]

% Note: t_a = 0 in the original solution of Tsai et al. (1994)
t_a = 0.0003;  % [m]
E_a = 3.1e9; % [Pa]
G_a = 1.1e9; % [Pa]

% Initial adhesive crack length
la0 = 0;  % [m]

%% Basic Calculations

% See notes for the FBD and load signs

% Total length
L   = L_1 + L_2;
% Initital geometry
l_a = L_1-(L_2-la0);
l_b = L_2-la0;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t + t_a)/(2*(l_a+l_b));

% Horizontal and vertical applied load component
Px = P*cos(alpha);
Py = P*sin(alpha);

% Bending stiffness - free adherent
D_11 = E*t^3/12;
% Bending stiffness - overlap region (lumped together)
D_00 = 2*E*(t^3/12+t*(t/2+t_a/2)^2) + E_a*t_a^3/12;

% x-vector spanning the free adherent length
x1 = 0:l_a/1000:l_a;
% x-vector spanning the overlap length
x0 = 0:l_b/1000:l_b;
% x-vector spanning the entire length
x = [x1 x0+l_a];

%% Roller-Roller
% Equivalent to Lai et al (1990), but self derrived

% Integration constants (see Maple for derrivations)

% Diff. eq. solution eigenvalues
lambda_1 = sqrt(P/D_11);
lambda_0 = sqrt(P/D_00);

A_1 = 0;
B_1 = -sqrt(D_11) * (cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * t + cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * t_a + (2 * alpha * l_a) + 0.2e1 * alpha * l_b - t - t_a) / (sqrt(D_11) * cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) + sqrt(D_00) * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b)) / 0.2e1;
A_0 = (-0.2e1 * sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * D_11 ^ (-0.1e1 / 0.2e1) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) + (-0.2e1 * l_a - 0.2e1 * l_b) * alpha + t + t_a) / cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) / 0.2e1;
B_0 = sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * D_11 ^ (-0.1e1 / 0.2e1);

R_A = 0;
M_A = 0;

% Vertical displacement of the free adherent
w1 = A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha + R_A/P) .* x1 + M_A/P;
% Vertical desiplacement of the overlap region
w0 = A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + alpha .* (l_a + x0) - (t + t_a)/2 + (l_a + x0) * R_A/P + M_A/P;

% Bending moment of the free adherent
M1 = P * (-A_1 * cosh(lambda_1 .* x1) - B_1 * sinh(lambda_1 .* x1));
% Bending moment of the overlap region
M0 = P * (-A_0 * cosh(lambda_0 .* x0) - B_0 * sinh(lambda_0 .* x0));

% Shear force in the free adherent
Q1 = -D_11*(A_1*P^(3/2)*sinh(sqrt(P)*x1/sqrt(D_11))/D_11^(3/2)+B_1*P^(3/2)*cosh(sqrt(P)*x1/sqrt(D_11))/D_11^(3/2));
% Shear force in the overlap region
Q0 = -D_00*(A_0*P^(3/2)*sinh(sqrt(P)*x0/sqrt(D_00))/D_00^(3/2)+B_0*P^(3/2)*cosh(sqrt(P)*x0/sqrt(D_00))/D_00^(3/2));

%% Clamped - Clamped

Method = 1;
if Method == 1
    % >> Geometric non linear solution & self derrived
    % >> Step 1: displacement and forces based on Lai et al (1990)
    
    % Diff. eq. solution eigenvalues
    lambda_1 = sqrt(P/D_11)
    lambda_0 = sqrt(P/D_00)

    % Integration constant
    a_mat = [0 0 0.1e1 / P 0 1 0 0 0; 0.1e1 / P 0 0 0 0 D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) 0 0; 0 0.1e1 / P * (l_b + l_a) 0.1e1 / P 0 0 0 cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b); 0 0.1e1 / P 0 0 0 0 D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b); 0.1e1 / P * l_a -0.1e1 / P * l_a 0 0 cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) -1 0; 0.1e1 / P -0.1e1 / P 0 0 D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) 0 -D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P); 1 -1 0 0 0 0 0 0; -l_b - l_a 0 -1 1 0 0 0 0]; 
    %a_mat = [0 0 0.1e1 / P 0 1 0 0 0; 0.1e1 / P 0 0 0 0 lambda_1 0 0; 0 0.1e1 / P * (l_b + l_a) 0.1e1 / P 0 0 0 cosh(lambda_0 * l_b) sinh(lambda_0 * l_b); 0 0.1e1 / P 0 0 0 0 lambda_0 * sinh(lambda_0 * l_b) lambda_0 * cosh(lambda_0 * l_b); 0.1e1 / P * l_a -0.1e1 / P * l_a 0 0 cosh(lambda_1 * l_a) sinh(lambda_1 * l_a) -1 0; 0.1e1 / P -0.1e1 / P 0 0 lambda_1 * sinh(lambda_1 * l_a) lambda_1 * cosh(lambda_1 * l_a) 0 -lambda_0; 1 -1 0 0 0 0 0 0; -l_b - l_a 0 -1 1 0 0 0 0];
    d_mat = [0 -alpha -(l_b + l_a) * alpha + t / 0.2e1 + t_a / 0.2e1 -alpha -t / 0.2e1 - t_a / 0.2e1 0 0 0]';
    b_mat = a_mat^(-1)*d_mat;
    
    % Isolote integration constants
    R_A = b_mat(1);
    R_B = b_mat(2);
    M_A = b_mat(3);
    M_B = b_mat(4);
    A_1 = b_mat(5);
    B_1 = b_mat(6);
    A_0 = b_mat(7);
    B_0 = b_mat(8);
  
    b_mat
    
    % Vertical displacement of the free adherent
    w11 = A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha + R_A/P) .* x1 + M_A/P;
    % Vertical desiplacement of the overlap region
    w00 = A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + alpha .* (l_a + x0) - (t + t_a)/2 + (l_a + x0) * R_A/P + M_A/P;
    
    % Bending moment in the free adherent
    M11 = P * (-A_1 * cosh(lambda_1 .* x1) - B_1 * sinh(lambda_1 .* x1));
    % Bending moment in the overlap region
    M00 = P * (-A_0 * cosh(lambda_0 .* x0) - B_0 * sinh(lambda_0 .* x0));
    
    % Shear force in the free adherent
    Q11 = -D_11*(A_1*P^(3/2)*sinh(sqrt(P) .* x1/sqrt(D_11))/D_11^(3/2)+B_1*P^(3/2)*cosh(sqrt(P) .* x1/sqrt(D_11))/D_11^(3/2));
    % Shear force in the overlap region
    Q00 = -D_00*(A_0*P^(3/2)*sinh(sqrt(P) .* x0/sqrt(D_00))/D_00^(3/2)+B_0*P^(3/2)*cosh(sqrt(P) .* x0/sqrt(D_00))/D_00^(3/2));
    
    % Step 2: Stresses based on F. Edde (1990)
    M_k = M11(1);
    Q_k = Q00(1);
    t   = t;
    E   = E;
    l_b = l_b;
    
    lambda = sqrt(8*G_a/(t*E*t_a));
    
    a_1 = sqrt(G_a)*(P*t+6*M_k)*sqrt(2)/(4*t^(3/2)*sqrt(E)*sqrt(t_a));
    a_2 = -(3*G_a*sqrt(t)*sqrt(E)*sqrt(t_a)*sqrt(2)*R_B+2*sinh(2*sqrt(2)*sqrt(G_a)*l_b/(sqrt(t)*sqrt(E)*sqrt(t_a)))*G_a^(3/2)*P*t+12*M_k*sinh(2*sqrt(2)*sqrt(G_a)*l_b/(sqrt(t)*sqrt(E)*sqrt(t_a)))*G_a^(3/2))*sqrt(2)/(8*G_a*cosh(2*sqrt(2)*sqrt(G_a)*l_b/(sqrt(t)*sqrt(E)*sqrt(t_a)))*t^(3/2)*sqrt(E)*sqrt(t_a));
    a_3 = (3*(R_B))/(4*t);
    
    Shear_a = a_1*sinh(lambda.*x0)+a_2*cosh(lambda.*x0)+a_3;
    
elseif Method == 2
    % >> Geometric Linear & original derrivation from F. Edde (1990)
    
    E = E;
    t = t;
    v = v;
    G = G_a;
    l_a = l_a;
    l_b = l_b;
    
    lambda = ((1+3*(1-v^2))*2*G_a/(E*t_a*t))^(1/2);
    
    a1 = 6*G*(-v^2+1)/(t_a*E*t^2);
    a2 = P*G/(t_a*E*t);
    a3 = 6*G*(-v^2+1)*t_a*E*t^2*lambda^2;
    a4 = t*l_b/lambda-(1/2)*t*l_b^2-t/lambda^2;
    a5 = (1/2)*l_b^2-l_a^2;
    a6 = (1/6)*t*l_b^3+t*l_b/lambda^2-t*l_b^2/(2*lambda)-t/lambda^3;
    a7 = (1/6)*l_b^3+(1/3)*l_a^3;
    
    % Solve for the 6 unknowns (see paper and Maple)
    a_mat = [lambda 0 0 -a1 * l_a -6 * G * (-v ^ 2 + 1) / t_a / E / t ^ 2 0; 0 1 0 a3 0 0; 1 / lambda 1 / lambda - 1 1 0 0 0; 0 0 -t - t_a -l_a - l_b -1 1; -t / lambda ^ 2 -a4 0 a5 -2 * l_a -l_b; t / lambda ^ 3 -a6 0 -a7 -2 * l_a ^ 2 l_b ^ 2 / 0.2e1];
    d_mat = [-P * G / t_a / E / t 0 0 0 0 0]';
    %b_mat = pinv(a_mat)*d_mat;
    b_mat = a_mat\d_mat;
    
    % Isolote integration constants
    B_1 = b_mat(1);
    B_2 = b_mat(2);
    P_1 = b_mat(3);
    Q   = b_mat(4);
    C   = b_mat(5);
    C_12 = b_mat(6);
    
    B_1 = -4*E*(l_a+(1/2)*l_b)*t_a*a2*((-l_a*l_b^2+(l_a*t*a3+(l_a*t_a-(1/2)*a4)*a3-l_a^2-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*lambda-l_a*l_b*a3*(t+t_a))*lambda^3*t^2/(4*E*(-l_a*l_b^2+(l_a*t*a3+(l_a*t_a-(1/2)*a4)*a3-l_a^2-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*(l_a+(1/2)*l_b)*t_a*t^2*lambda^5-4*E*(l_a+(1/2)*l_b)*t_a*l_b*l_a*a3*(t+t_a)*t^2*lambda^4+4*l_b*(((1/2)*E*l_a^2*t^2*t_a*a1-3*G*(v-1)*(v+1)*(a3*a4+a5)*(1/2))*l_b+E*l_a^3*t^2*t_a*a1+3*G*(v-1)*(v+1)*(-a3*a6+a7))*(t+t_a)*lambda^3-(4*(G*(-3*v^2*(1/2)+3/2)*l_b^3+(-(1/4)*E*l_a*t^2*t_a*a1+3*G*a3*(v-1)*(v+1)*t*(1/2)-3*G*(v-1)*(v+1)*(-a3*t_a+l_a)*(1/2))*l_b^2+E*l_a^3*t^2*t_a*a1+3*G*(v-1)*(v+1)*(-a3*a6+a7)))*t*lambda^2-(4*(-3*G*(v-1)*(v+1)*(a3*t+a3*t_a-2)*l_b^2*(1/2)+((1/2)*E*l_a*t^2*t_a*a1+(-3*v^2+3)*G*a3*t+3*G*(v-1)*(v+1)*(-a3*t_a+l_a))*l_b+E*l_a^2*t^2*t_a*a1-3*G*(v-1)*(v+1)*(a3*a4+a5)))*t*lambda-12*l_b*t*G*a3*(v-1)*(v+1)*(t+t_a));
    B_2 = -4*E*(l_a+(1/2)*l_b)*t_a*a2*a3*(l_a*l_b*(t+t_a)*lambda^2-(l_a-(1/2)*l_b)*t*lambda-t)*lambda*t^2/(4*E*(l_b*l_a*t*a3-l_a*l_b^2+(a3*l_a*t_a-l_a^2-(1/2)*a3*a4-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*(l_a+(1/2)*l_b)*t_a*t^2*lambda^5-4*E*(l_a+(1/2)*l_b)*t_a*l_b*l_a*a3*(t+t_a)*t^2*lambda^4+4*l_b*(E*l_a^2*t_a*a1*(l_a+(1/2)*l_b)*t^2+(3*(v+1))*(v-1)*G*((-(1/2)*a3*a4-(1/2)*a5)*l_b-a3*a6+a7))*(t+t_a)*lambda^3-(4*((l_a-(1/2)*l_b)*E*(l_a+(1/2)*l_b)*t_a*l_a*a1*t^2+3*l_b^2*G*a3*(v-1)*(v+1)*t*(1/2)+(3*(v+1))*(v-1)*G*(-(1/2)*l_b^3+((1/2)*a3*t_a-(1/2)*l_a)*l_b^2-a3*a6+a7)))*t*lambda^2-(4*(E*l_a*t_a*a1*(l_a+(1/2)*l_b)*t^2-3*l_b*G*a3*(v-1)*(v+1)*(l_b+2)*t*(1/2)+(3*(v+1))*((-(1/2)*a3*t_a+1)*l_b^2+(-a3*t_a+l_a)*l_b-a3*a4-a5)*(v-1)*G))*t*lambda-12*l_b*t*G*a3*(v-1)*(v+1)*(t+t_a));
    C   = -2*E*t_a*(l_b*((-(1/2)*a3*a4-(1/2)*a5)*l_b-a3*a6+a7)*(t+t_a)*lambda^3-((1/2)*l_b^2*t*a3-(1/2)*l_b^3+((1/2)*a3*t_a-(1/2)*l_a)*l_b^2-a3*a6+a7)*t*lambda^2-(-(1/2)*l_b*a3*(l_b+2)*t+(-(1/2)*a3*t_a+1)*l_b^2+(-a3*t_a+l_a)*l_b-a3*a4-a5)*t*lambda-l_b*t*a3*(t+t_a))*a2*t^2/(4*E*(l_b*l_a*t*a3-l_a*l_b^2+((l_a*t_a-(1/2)*a4)*a3-l_a^2-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*(l_a+(1/2)*l_b)*t_a*t^2*lambda^5-4*E*(l_a+(1/2)*l_b)*t_a*l_b*l_a*a3*(t+t_a)*t^2*lambda^4+4*l_b*(E*l_a^2*t_a*a1*(l_a+(1/2)*l_b)*t^2+(3*(v+1))*(v-1)*G*((-(1/2)*a3*a4-(1/2)*a5)*l_b-a3*a6+a7))*(t+t_a)*lambda^3-(4*((l_a-(1/2)*l_b)*E*(l_a+(1/2)*l_b)*t_a*l_a*a1*t^2+3*l_b^2*G*a3*(v-1)*(v+1)*t*(1/2)+(3*(v+1))*(v-1)*G*(-(1/2)*l_b^3+((1/2)*a3*t_a-(1/2)*l_a)*l_b^2-a3*a6+a7)))*t*lambda^2-(4*(E*l_a*t_a*a1*(l_a+(1/2)*l_b)*t^2-3*l_b*G*a3*(v-1)*(v+1)*(l_b+2)*t*(1/2)+(3*(v+1))*((-(1/2)*a3*t_a+1)*l_b^2+(-a3*t_a+l_a)*l_b-a3*a4-a5)*(v-1)*G))*t*lambda-12*l_b*t*G*a3*(v-1)*(v+1)*(t+t_a));
    C_12 = 4*E*(l_a*((a4*l_a-a6)*a3+l_a*a5+a7)*(t+t_a)*lambda^3+(1/2)*(2*l_a^2*t*a3-2*l_b*l_a^2+(2*l_a^2*t_a-a6)*a3-2*l_a^3+a7)*t*lambda^2-(l_a*a3*(l_a-1)*t+l_a*l_b+((l_a^2-l_a)*t_a+(1/2)*a4)*a3+l_a^2+(1/2)*a5)*t*lambda-l_a*t*a3*(t+t_a))*t_a*a2*t^2/(4*E*(l_b*l_a*t*a3-l_a*l_b^2+((l_a*t_a-(1/2)*a4)*a3-l_a^2-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*(l_a+(1/2)*l_b)*t_a*t^2*lambda^5-4*E*(l_a+(1/2)*l_b)*t_a*l_b*l_a*a3*(t+t_a)*t^2*lambda^4+4*l_b*(E*l_a^2*t_a*a1*(l_a+(1/2)*l_b)*t^2+(3*(v+1))*(v-1)*G*((-(1/2)*a3*a4-(1/2)*a5)*l_b-a3*a6+a7))*(t+t_a)*lambda^3-(4*((l_a-(1/2)*l_b)*E*(l_a+(1/2)*l_b)*t_a*l_a*a1*t^2+3*l_b^2*G*a3*(v-1)*(v+1)*t*(1/2)+(3*(v+1))*(v-1)*G*(-(1/2)*l_b^3+((1/2)*a3*t_a-(1/2)*l_a)*l_b^2-a3*a6+a7)))*t*lambda^2-(4*(E*l_a*t_a*a1*(l_a+(1/2)*l_b)*t^2-3*l_b*G*a3*(v-1)*(v+1)*(l_b+2)*t*(1/2)+(3*(v+1))*((-(1/2)*a3*t_a+1)*l_b^2+(-a3*t_a+l_a)*l_b-a3*a4-a5)*(v-1)*G))*t*lambda-12*l_b*t*G*a3*(v-1)*(v+1)*(t+t_a));
    P_1 = 4*E*(l_a+(1/2)*l_b)*t_a*a2*((-l_a*l_b^2+(-l_a^2-(1/2)*a3*a4-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*lambda^3+t*a3*(l_a-(1/2)*l_b)*lambda^2-(l_a-(1/2)*l_b-1)*a3*t*lambda-a3*t)*t^2/(4*E*(l_b*l_a*t*a3-l_a*l_b^2+((l_a*t_a-(1/2)*a4)*a3-l_a^2-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*(l_a+(1/2)*l_b)*t_a*t^2*lambda^5-4*E*(l_a+(1/2)*l_b)*t_a*l_b*l_a*a3*(t+t_a)*t^2*lambda^4+4*l_b*(E*l_a^2*t_a*a1*(l_a+(1/2)*l_b)*t^2+(3*(v+1))*(v-1)*G*((-(1/2)*a3*a4-(1/2)*a5)*l_b-a3*a6+a7))*(t+t_a)*lambda^3-(4*((l_a-(1/2)*l_b)*E*(l_a+(1/2)*l_b)*t_a*l_a*a1*t^2+3*l_b^2*G*a3*(v-1)*(v+1)*t*(1/2)+(3*(v+1))*(v-1)*G*(-(1/2)*l_b^3+((1/2)*a3*t_a-(1/2)*l_a)*l_b^2-a3*a6+a7)))*t*lambda^2-(4*(E*l_a*t_a*a1*(l_a+(1/2)*l_b)*t^2-3*l_b*G*a3*(v-1)*(v+1)*(l_b+2)*t*(1/2)+(3*(v+1))*((-(1/2)*a3*t_a+1)*l_b^2+(-a3*t_a+l_a)*l_b-a3*a4-a5)*(v-1)*G))*t*lambda-12*l_b*t*G*a3*(v-1)*(v+1)*(t+t_a));
    Q   = 4*E*(l_a+(1/2)*l_b)*t_a*a2*(l_a*l_b*(t+t_a)*lambda^2-(l_a-(1/2)*l_b)*t*lambda-t)*lambda*t^2/(4*E*(l_b*l_a*t*a3-l_a*l_b^2+(a3*l_a*t_a-l_a^2-(1/2)*a3*a4-(1/2)*a5)*l_b+(a4*l_a-a6)*a3+l_a*a5+a7)*(l_a+(1/2)*l_b)*t_a*t^2*lambda^5-4*E*(l_a+(1/2)*l_b)*t_a*l_b*l_a*a3*(t+t_a)*t^2*lambda^4+4*l_b*(E*l_a^2*t_a*a1*(l_a+(1/2)*l_b)*t^2+(3*(v+1))*(v-1)*G*((-(1/2)*a3*a4-(1/2)*a5)*l_b-a3*a6+a7))*(t+t_a)*lambda^3-(4*((l_a-(1/2)*l_b)*E*(l_a+(1/2)*l_b)*t_a*l_a*a1*t^2+3*l_b^2*G*a3*(v-1)*(v+1)*t*(1/2)+(3*(v+1))*(v-1)*G*(-(1/2)*l_b^3+((1/2)*a3*t_a-(1/2)*l_a)*l_b^2-a3*a6+a7)))*t*lambda^2-(4*(E*l_a*t_a*a1*(l_a+(1/2)*l_b)*t^2-3*l_b*G*a3*(v-1)*(v+1)*(l_b+2)*t*(1/2)+(3*(v+1))*((-(1/2)*a3*t_a+1)*l_b^2+(-a3*t_a+l_a)*l_b-a3*a4-a5)*(v-1)*G))*t*lambda-12*l_b*t*G*a3*(v-1)*(v+1)*(t+t_a));
    
    B_1
    B_2
    P_1
    Q
    C
    C_12
    
    A_1 = -B_2;
    A_2 = B_2;
    
    Shear_a = A_1*sinh(lambda*x0)+A_2*(cosh(lambda*x0)-1);
    
end

figure(1)
hold on
% plot(x, -[w1 w0], 'g')
% plot(x(length(w1)+1:end), -w0, 'r')
plot(x, -[w11 w00], 'b')
plot(x(length(w11)+1:end), -w00, 'c')
hold off

figure(2)
hold on
% plot(x, -[M1 M0], 'g')
% plot(x(length(M1)+1:end), -M0, 'r')
plot(x, [M11 M00], 'b')
plot(x(length(M11)+1:end), M00, 'c')
hold off
% 
% figure(3)
% hold on
% plot(x, -[Q1 Q0], 'g')
% plot(x(length(Q1)+1:end), -Q0, 'r')
% plot(x, -[Q11 Q00], 'b')
% plot(x(length(Q11)+1:end), -Q00, 'c')
% hold off
% 
% figure(4)
% hold on
% plot(x0, Shear_a, 'g')
% hold off


