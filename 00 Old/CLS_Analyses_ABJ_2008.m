function [w1, w0, M1, M0, x_load, shear_a, peel_a, x_stress] = CLS_Analyses_ABJ_2008(Load, Strap, Lap, Adhesive, Conditions)
% Author    :   Modeling of Adhesive Bonded Joints > Luo and Tong
% Method    :   Goland and Reissner
% Remarks   :   Stresses added as they are not in the paper
%           :   Integration constants solutions from the book
%           :   Compared and confirmed with solution from the self derrived
%               solution of Lai et al. (1994).

% Inputs
P   = Load.P;
R_A = Load.R_A;
M_A = Load.M_A;

E_1 = Strap.E;
t_1 = Strap.t;
v_1 = Strap.v;
L_1 = Strap.L;

E_2 = Lap.E;
t_2 = Lap.t;
v_2 = Lap.v;
L_2 = Lap.L;

t_a = Adhesive.t;
E_a = Adhesive.E;
G_a = Adhesive.G;
%t_a = 0;

% Define geometry
L   = L_1;
l_a = (L - L_2) + Adhesive.la0;
l_b = L_2 - Adhesive.la0;

l = l_a
c = l_b

alpha = (t_2 + t_a)/(2*(l+c));

% Note: F in this solution is the horizontal component of the load P (P =
% applied along the neutral line axis)
F = P * (l+c)/ sqrt((l+c)^2+(t_a/2+t_2/2)^2)


% FOR NOW. THIS NEEDS WORK FOR GENERIC CROSS-SECTIONS
D_11 = E_1 * 1/12 * t_1^3;
D_00 = 2 * E_2 * (1/12 * t_1^3 + t_1 * (t_1/2)^2); % >>>> FOR NOW SYMMETRICA
A_11 = E_1 * t_1;
A_22 = E_1 * (t_1 + t_2);

F = (8/c)^2*D_11

beta_k = sqrt(F/D_11);
beta_0 = sqrt(F/D_00);

beta_k*c

% Y instead of A for the integration constants
Y_1 = -1 * (beta_0*(t_2 + t_a) * cosh(beta_0*c)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
B_1 = -1 * (beta_k*(t_2 + t_a) * cosh(beta_k*l)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
Y_0 = 0;
B_0 = 0;

x1 = 0:l_a/100:l_a;
x0 = -l_b:l_b/100:0;

x_load = [x1 x0+l_b+l_a];

% Vertical displacement
w1 = Y_1 * sinh(beta_k * x1) + alpha * x1;
w0 = B_1 * sinh(beta_0 * x0) + alpha * x0;

% Bending moment
M1 = -D_11 * Y_1 * beta_k^2 * sinh(beta_k .* x1);
M0 = -D_00 * B_1 * beta_0^2 * sinh(beta_0 .* x0);

% Adhesive shear stress
k = 1/(1+(beta_k/beta_0)*tanh(beta_0*c)*coth(beta_k*l))
k = 1/(1+sqrt(D_00/D_11)*tanh(sqrt(D_11/D_00)*beta_k*c)*coth(beta_k*l))

% % Find the integration constants
% alpha_k = A_11 * t_1^2 / (4 * D_11);
% alpha_a = (1 + alpha_k) / 4;
% 
% M_k = M1(end);
% 
% V_k = 0;
% 
% beta_s = sqrt(2)/2 * (2*E_a / (D_11*t_a))^(1/4);
% beta_t = sqrt(8 * G_a / (A_11*t_a));
% beta_a = sqrt(alpha_a * beta_t^2);
% 
% beta_s1 = M_k*(sinh(beta_s*c*cos(beta_s*c))+cosh(beta_s*c)*sin(beta_0*c)) + V_k/beta_s*sinh(beta_s*c)*sin(beta_s*c);
% beta_s4 = M_k*(sinh(beta_s*c*cos(beta_s*c))-cosh(beta_s*c)*sin(beta_0*c)) + V_k/beta_s*cosh(beta_s*c)*cos(beta_s*c);
% 
% % x-steps for adhesive stresses
% 
% % Compute the stresses
% Shear_a = beta_t*(F_t1+6*M_k)*cosh(beta_t.x)/(8*t_1*sinh(beta_t*c) + 3*(F_t1-2*M_k)/(8*t_1*c));
% Peel_a = 2*beta_s^2*(B_s1*sinh(beta_s*x)*sin(beta_s*x)+B_s4*cosh(beta_s*x)*cos(beta_s*x)) / (sinh(2*beta_s*c)+sin(2*beta_s));


end