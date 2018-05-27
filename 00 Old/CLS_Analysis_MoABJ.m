clear all; close all; clc

% Load
S_max = 200e8; % [N/m^2]
R_load = 0.1; % [-]

% Total number of cycles
N = 1e9; % [-]

% Adherent
E   = 70e9; % [Pa]
t   = 0.0016; % [m]
v   = 0.34; % [-]
% Strap (upper, long adherent) length
L_1 = 32*0.0016*(1+1.25); % [m]
% lap (lower, short adherent) length
L_2 = 32*0.0016; % [m]

% Adhesive
t_a = 0.078*0.0016; % [m]
E_a = 0.04*70e9; % [Pa]
G_a = 1.1e9; % [Pa]
v_a = 0.4; % [-]

% Initial crack length
b = 0; % [m]

%% Geometry

% See notes for the FBD and load signs

% Total length
L   = L_1;
% Overlap length
l_A = (L - L_2) +b;
% Free adherent length
l_B = L_2 - b;

% Re-Define to match the definitions used in the book
l = l_A;
c = l_B;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t + t_a)/(2*(l+c));

% Adherent stiffnesses
% Free adherent
D_11 = E*t^3/12;
% Overlap region (lumped together)
D_00 = 2*E*(t^3/12+t*(t/2+t_a/2)^2) + E_a*t_a^3/12;

% Axial stiffness - free adherent
A_11 = E * t + E_a*t_a;
% Axial stiffness - overlap region (lumped together)
% Note: adhesive ignored
A_22 = E * (t + t) + E_a*t_a;

%% Applied load

% Note: F in this solution is the horizontal component of the load P (P =
% applied along the neutral line axis)
P=(8/l_B)^2*D_11;
F = P * (l+c)/ sqrt((l+c)^2+(t_a/2+t/2)^2);

% To compare with plots in the book
%F = (8/c)^2*D_11;

%% Axis preperation

% x-vector spanning the free adherent length 
x1 = 0:l_A/1000:l_A;
% x-vector spanning the overlap length
x0 = -l_B:l_B/1000:0;
% x-vector spanning the entire length
x = [x1 x0+l_B+l_A];

%% Step 1 - Displacements and forces

% Roller-Roller BC
beta_k = sqrt(F/D_11);
beta_0 = sqrt(F/D_00);

% Y instead of A for the integration constants
Y_1 = -1 * (beta_0*(t + t_a) * cosh(beta_0*c)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
B_1 = -1 * (beta_k*(t + t_a) * cosh(beta_k*l)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
Y_0 = 0;
B_0 = 0;

% Vertical displacement of the free adherent
w1 = Y_1 * sinh(beta_k * x1) + alpha * x1;
% Vertical desiplacement of the overlap region
w0 = B_1 * sinh(beta_0 * x0) + alpha * x0;

% Bending moment of the free adherent
M1 = -D_11 * Y_1 * beta_k^2 * sinh(beta_k .* x1);
% Bending moment of the overlap region
M0 = -D_00 * B_1 * beta_0^2 * sinh(beta_0 .* x0);

Q1 = -D_11 * Y_1 * beta_k^3 * cosh(beta_k*x1);
Q0 = -D_00 * B_1 * beta_0^3 * cosh(beta_0*x0);
        
% Edge Moment Factor
%k = 1/(1+(beta_k/beta_0)*tanh(beta_0*c)*coth(beta_k*l));
%k = 1/(1+sqrt(D_00/D_11)*tanh(sqrt(D_11/D_00)*beta_k*c)*coth(beta_k*l));

%% Step 2 - Stresses

% Find the integration constants
alpha_k = A_11 * t^2 / (4 * D_11);
alpha_a = (1 + alpha_k) / 4;
 
M_k = M1(end)
 
V_k = 0;
V_k = Q1(end)
 
beta_s = sqrt(2)/2 * (2*E_a / (D_11*t_a))^(1/4);
beta_t = sqrt(8 * G_a / (A_11*t_a));
beta_a = sqrt(alpha_a * beta_t^2);

B_s1 = M_k*(sinh(beta_s*c*cos(beta_s*c))+cosh(beta_s*c)*sin(beta_0*c)) + V_k/beta_s*sinh(beta_s*c)*sin(beta_s*c);
B_s4 = M_k*(sinh(beta_s*c*cos(beta_s*c))-cosh(beta_s*c)*sin(beta_0*c)) + V_k/beta_s*cosh(beta_s*c)*cos(beta_s*c);

% Adhesive shear stress
shear_a = beta_t*(F*t+6*M_k)*cosh(beta_t*x0)/(8*t*sinh(beta_t*c) + 3*(F*t-2*M_k)/(8*t*c));
% Adheisve peel stress
peel_a = 2*beta_s^2*(B_s1*sinh(beta_s*x0).*sin(beta_s*x0)+B_s4*cosh(beta_s*x0).*cos(beta_s*x0)) / (sinh(2*beta_s*c)+sin(2*beta_s*c));

figure(1)
plot(x0/l_B,-w0/t)

figure(2)
plot(x0(1:end/2)/l_B,shear_a(1:end/2)/10^6)

figure(3)
plot(x0(1:end/2)/l_B,peel_a(1:end/2)/10^6)




% Clamped-Clamped BC
%A_1 = beta_0*(2*cosh(beta_0*c)*sinh(beta_0*c)*beta_0*t+2*cosh(beta_0*c)*sinh(beta_0*c)*beta_0*t_a+cosh(beta_0*c)*sinh(beta_k*l)*beta_k*t+cosh(beta_0*c)*sinh(beta_k*l)*beta_k*t_a-cosh(beta_k*l)*sinh(beta_0*c)*beta_0*t-cosh(beta_k*l)*sinh(beta_0*c)*beta_0*t_a-beta_0*sinh(beta_0*c)*t-beta_0*sinh(beta_0*c)*t_a-sinh(beta_k*l)*beta_k*t-sinh(beta_k*l)*beta_k*t_a)/(2*(cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*c+cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*l+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*c+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*l-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*c+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_0*c)^2*beta_0*beta_k-2*cosh(beta_0*c)*cosh(beta_k*l)*beta_0*beta_k+cosh(beta_k*l)^2*beta_0*beta_k-sinh(beta_0*c)^2*beta_0*beta_k-sinh(beta_0*c)*sinh(beta_k*l)*beta_0^2-sinh(beta_0*c)*sinh(beta_k*l)*beta_k^2-sinh(beta_k*l)^2*beta_0*beta_k));
%A_2 = (-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c*t-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c*t_a-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l*t-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l*t_a-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*c*t-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*c*t_a-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*l*t-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*l*t_a+cosh(beta_0*c)^2*beta_k*t+cosh(beta_0*c)^2*beta_k*t_a-cosh(beta_0*c)*cosh(beta_k*l)*beta_k*t-cosh(beta_0*c)*cosh(beta_k*l)*beta_k*t_a+sinh(beta_0*c)^2*beta_k*t+sinh(beta_0*c)^2*beta_k*t_a+sinh(beta_0*c)*sinh(beta_k*l)*beta_0*t+sinh(beta_0*c)*sinh(beta_k*l)*beta_0*t_a+sinh(beta_0*c)*beta_k*c*t+sinh(beta_0*c)*beta_k*c*t_a+sinh(beta_0*c)*beta_k*l*t+sinh(beta_0*c)*beta_k*l*t_a-cosh(beta_0*c)*beta_k*t-cosh(beta_0*c)*beta_k*t_a+cosh(beta_k*l)*beta_k*t+cosh(beta_k*l)*beta_k*t_a)*beta_0/(2*(cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*c+cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*l+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*c+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*l-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*c+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_0*c)^2*beta_0*beta_k-2*cosh(beta_0*c)*cosh(beta_k*l)*beta_0*beta_k+cosh(beta_k*l)^2*beta_0*beta_k-sinh(beta_0*c)^2*beta_0*beta_k-sinh(beta_0*c)*sinh(beta_k*l)*beta_0^2-sinh(beta_0*c)*sinh(beta_k*l)*beta_k^2-sinh(beta_k*l)^2*beta_0*beta_k));

% Avoids the MATLAB rounding error
% a11 = beta_k*beta_0*(l+c)*cosh(beta_0*c)-beta_k*sinh(beta_0*c)-beta_0*sinh(beta_k*l);
% a12 = beta_0*(cosh(beta_0*c)-cosh(beta_k*l));
% c1  = (t+t_a)/2*beta_0*(cosh(beta_0*c)-1);
% 
% a21 = beta_k*(cosh(beta_0*c)-cosh(beta_k*l))-beta_k*(l+c)*sinh(beta_0*c);
% a22 = -1*(beta_k*sinh(beta_k*l)+beta_0*sinh(beta_0*c));
% c2  = (t+t_a)/2*beta_0*sinh(beta_0*c);
% 
% Values =  [a11 a12 ; a21 a22]^(-1) * [c1 ; c2];
% 
% A_1 = Values(1);
% A_2 = Values(2);
% 
% B_1 = beta_k/beta_0*A_1;
% B_2 = A_2+A_1*beta_k*(l+c)+(t+t_a)/2;
% 
% M_A = -A_2*F;
% R_A = -beta_k*A_1*F;
% M_B = -B_2*F;
% R_B = beta_0*B_1*F;
% 
% M_B - M_A + F*(t + t_a)/2 - R_A*(l+c)
% 
% w1 = A_1 * sinh(beta_k*x1) + A_2*cosh(beta_k*x1) + R_A/F*x1 + M_A/F;
% w0 = B_1 * sinh(beta_0*x0) + B_2*cosh(beta_0*x0) - R_B/F*x0 + M_B/F;

%w11 = A_1*sqrt(F)*sinh(sqrt(F)*x_1/sqrt(D_11))/sqrt(D_11)+B_1*sqrt(F)*cosh(sqrt(F)*x_1/sqrt(D_11))/sqrt(D_11)-R_A/F;
%w22 = A_0*sqrt(F)*sinh(sqrt(F)*x_0/sqrt(D_2))/sqrt(D_2)+B_0*sqrt(F)*cosh(sqrt(F)*x_0/sqrt(D_2))/sqrt(D_2)+R_A/F;

% figure(1)
% hold on
% plot(x, -[w1 w0]/t, 'g')
% plot(x(length(w1)+1:end), -w0/t, 'r')
% hold off

% w11 = A_1*(sinh(beta_k*x1)-beta_k*x1)+A_2*(cosh(beta_k*x1)-1);
% w00 = A_1*(beta_k/beta_0*(sinh(beta_0*x0)-beta_0*x0)+ ...
%             beta_k*(l+c)*(cosh(beta_0*x0)-1)) ...
%             +A_2*(cosh(beta_0*x0)-1)+(t+t_a)/2*(cosh(beta_0*x0)-1);