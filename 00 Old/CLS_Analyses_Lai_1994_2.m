function [w1, w0, M1, M0, x] = CLS_Analyses_Lai_1994_2(Load, Strap, Lap, Adhesive, Conditions)
% Author    :   Lai et al. (1994)
% Method    :   Goland and Reissner
% Remarks   :   Stresses added as they are not in the paper
%           :   Integration constants self derived
%           :   Compared and confirmed with solution from Modeling of
%               Adhesive Bonded Joints

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
%t_a = 0;

R_A = 0;
M_A = 0;

% Define geometry
L   = L_1;
l_a = (L - L_2) + Adhesive.la0;
l_b = L_2 - Adhesive.la0;

D_1 = E_1 * 1/12 * t_1^3;
D_0 = 2 * E_2 * (1/12 * t_1^3 + t_1 * (t_1/2)^2); % >>>> FOR NOW SYMMETRICA

lambda_1 = sqrt(P/D_1);
lambda_0 = sqrt(P/D_0);

alpha = (t_2 + t_a)/(2*(l_a+l_b));

A_1 = 0;
B_1 = -sqrt(D_1) * (cosh(D_0 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * t_2 + cosh(D_0 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * t_a + (2 * alpha * l_a) + 0.2e1 * alpha * l_b - t_2 - t_a) / (sqrt(D_1) * cosh(D_0 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * sinh(D_1 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) + sqrt(D_0) * cosh(D_1 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * sinh(D_0 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b)) / 0.2e1;
A_0 = (-0.2e1 * sqrt(D_0) * B_1 * cosh(D_1 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * D_1 ^ (-0.1e1 / 0.2e1) * sinh(D_0 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) + (-0.2e1 * l_a - 0.2e1 * l_b) * alpha + t_2 + t_a) / cosh(D_0 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) / 0.2e1;
B_0 = sqrt(D_0) * B_1 * cosh(D_1 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * D_1 ^ (-0.1e1 / 0.2e1);

x1 = 0:l_a/100:l_a;
x0 = 0:l_b/100:l_b;

x = [x1 x0+l_a];

w1 = A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha + R_A/P) .* x1 + M_A/P;
w0 = A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + alpha .* (l_a + x0) - (t_1 + t_a)/2 + (l_a + x0) * R_A/P + M_A/P;

M1 = P * (-A_1 * cosh(lambda_1 .* x1) - B_1 * sinh(lambda_1 .* x1));
M0 = P * (-A_0 * cosh(lambda_0 .* x0) - B_0 * sinh(lambda_0 .* x0));

end