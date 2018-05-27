function [w1, w0, M1, M0, x] = CLS_Analysis_Lai_et_al_adjusted(Load, Strap, Lap, Adhesive, Conditions)
% Paper     > THE CRACKED LAP SHEAR SPECIMEN REVISITED - A CLOSED FORM
%           > SOLUTION
% Author    > Lai et al. (1994)
% Method    > Goland and Reissner
% Includes  > Normalized solutions for (1) overlap loads and
%           displacements, (2) SERR, (3) Mode Ratio and (4) Fracture
%           Efficieny
%           > Solutions for the BCs (1) C-C, (2) R-C, (3) R-R
%
% Remarks   > Adhesive stresses added as they are not included
%           > Adherent stresses added as they are not included
%           > Integration constants self derived; errors found in their
%           paper
%           > Compared and confirmed the custom derrivation with solution 
%           from 'Modeling of Adhesively Bonded Joints'
%% Load Input Data
P   = Load.P;

E_1 = Strap.E;
t_1 = Strap.t;
v_1 = Strap.v;
L_1 = Strap.L;

E_2 = Lap.E;
t_2 = Lap.t;
v_2 = Lap.v;
L_2 = Lap.L;

% Note: t_a = 0 in the original solution of Tsai et al. (1994)
t_a = Adhesive.t;
E_a = Adhesive.E;
G_a = Adhesive.G;

%% Geometry

% See notes for the FBD and load signs

% Total length
L   = L_1;
% Overlap length
l_a = (L - L_2) + Adhesive.la0;
% Free adherent length
l_b = L_2 - Adhesive.la0;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t_2 + t_a)/(2*(l_a+l_b));

% Bending stiffness - free adherent
D_11 = E_1 * 1/12 * t_1^3;
% Bending stiffness - overlap region (lumped together)
% To Do: include assymetric cross-sections
% To Do: include t_a
D_00 = 2 * E_2 * (1/12 * t_1^3 + t_1 * (t_1/2)^2);

%% Solution constants

% Diff. eq. solution eigenvalues
lambda_1 = sqrt(P/D_11);
lambda_0 = sqrt(P/D_00);

%% Axis preperation

% x-vector spanning the free adherent length 
x1 = 0:l_a/1000:l_a;
% x-vector spanning the overlap length
x0 = 0:l_b/1000:l_b;

% x-vector spanning the entire length
x = [x1 x0+l_a];

%% Solutions

% >>>>>>>>>>>>>>>>>>>>>>>>>>>> Roller-Roller <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Integration constants
A_1 = 0;
B_1 = -sqrt(D_11) * (cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * t_2 + cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * t_a + (2 * alpha * l_a) + 0.2e1 * alpha * l_b - t_2 - t_a) / (sqrt(D_11) * cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) * sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) + sqrt(D_00) * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b)) / 0.2e1;
A_0 = (-0.2e1 * sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * D_11 ^ (-0.1e1 / 0.2e1) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) + (-0.2e1 * l_a - 0.2e1 * l_b) * alpha + t_2 + t_a) / cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) / 0.2e1;
B_0 = sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) * D_11 ^ (-0.1e1 / 0.2e1);

R_A = 0;
M_A = 0;

% Vertical displacement of the free adherent
w1 = A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha + R_A/P) .* x1 + M_A/P;
% Vertical desiplacement of the overlap region
w0 = A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + alpha .* (l_a + x0) - (t_1 + t_a)/2 + (l_a + x0) * R_A/P + M_A/P;

% Bending moment of the free adherent
M1 = P * (-A_1 * cosh(lambda_1 .* x1) - B_1 * sinh(lambda_1 .* x1));
% Bending moment of the overlap region
M0 = P * (-A_0 * cosh(lambda_0 .* x0) - B_0 * sinh(lambda_0 .* x0));

% Shear force in the free adherent
Q1 = -D_11*(A_1*P^(3/2)*sinh(sqrt(P)*x1/sqrt(D_11))/D_11^(3/2)+B_1*P^(3/2)*cosh(sqrt(P)*x1/sqrt(D_11))/D_11^(3/2));
% Shear force in the overlap region
Q0 = -D_00*(A_0*P^(3/2)*sinh(sqrt(P)*x0/sqrt(D_00))/D_00^(3/2)+B_0*P^(3/2)*cosh(sqrt(P)*x0/sqrt(D_00))/D_00^(3/2));
% >>>>>>>>>>>>>>>>>>>>>>>> Clamped-Clamped <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
% Integration constant
a_mat = [0 0 0.1e1 / P 0 1 0 0 0; 0.1e1 / P 0 0 0 0 D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) 0 0; 0 0.1e1 / P * (l_b + l_a) 0.1e1 / P 0 0 0 cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b); 0 0.1e1 / P 0 0 0 0 D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b) D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_b); 0.1e1 / P * l_a -0.1e1 / P * l_a 0 0 cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) -1 0; 0.1e1 / P -0.1e1 / P 0 0 D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_a) 0 -D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P); 1 -1 0 0 0 0 0 0; -l_b - l_a 0 -1 1 0 0 0 0];
d_mat = [0 -alpha -(l_b + l_a) * alpha + t_2 / 0.2e1 + t_a / 0.2e1 -alpha -t_2 / 0.2e1 - t_a / 0.2e1 0 0 0]';
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

% Vertical displacement of the free adherent
w11 = A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha + R_A/P) .* x1 + M_A/P;
% Vertical desiplacement of the overlap region
w00 = A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + alpha .* (l_a + x0) - (t_1 + t_a)/2 + (l_a + x0) * R_A/P + M_A/P;

% Bending moment in the free adherent
M11 = P * (-A_1 * cosh(lambda_1 .* x1) - B_1 * sinh(lambda_1 .* x1));
% Bending moment in the overlap region
M00 = P * (-A_0 * cosh(lambda_0 .* x0) - B_0 * sinh(lambda_0 .* x0));

% Shear force in the free adherent
Q11 = -D_11*(A_1*P^(3/2)*sinh(sqrt(P) .* x1/sqrt(D_11))/D_11^(3/2)+B_1*P^(3/2)*cosh(sqrt(P) .* x1/sqrt(D_11))/D_11^(3/2));
% Shear force in the overlap region
Q00 = -D_00*(A_0*P^(3/2)*sinh(sqrt(P) .* x0/sqrt(D_00))/D_00^(3/2)+B_0*P^(3/2)*cosh(sqrt(P) .* x0/sqrt(D_00))/D_00^(3/2));

% Stresses
M_k = M11(1);
Q_k = Q00(1);
t   = t_1;
E   = E_1;
l_B = l_b;

lambda = sqrt(8*G_a/(t*E*t_a));

a_1 = sqrt(G_a)*(P*t+6*M_k)*sqrt(2)/(4*t^(3/2)*sqrt(E)*sqrt(t_a));
a_2 = -(3*G_a*sqrt(t)*sqrt(E)*sqrt(t_a)*sqrt(2)*R_B+2*sinh(2*sqrt(2)*sqrt(G_a)*l_B/(sqrt(t)*sqrt(E)*sqrt(t_a)))*G_a^(3/2)*P*t+12*M_k*sinh(2*sqrt(2)*sqrt(G_a)*l_B/(sqrt(t)*sqrt(E)*sqrt(t_a)))*G_a^(3/2))*sqrt(2)/(8*G_a*cosh(2*sqrt(2)*sqrt(G_a)*l_B/(sqrt(t)*sqrt(E)*sqrt(t_a)))*t^(3/2)*sqrt(E)*sqrt(t_a));
a_3 = (3*(R_B))/(4*t);

Shear_a = a_1*sinh(lambda.*x0)+a_2*cosh(lambda.*x0)+a_3;

figure(1)
hold on
plot(x, -[w1 w0], 'g')
plot(x(length(w1)+1:end), -w0, 'r')
plot(x, -[w11 w00], 'b')
plot(x(length(w11)+1:end), -w00, 'c')
hold off

% figure(2)
% hold on
% plot(x, -[M1 M0], 'g')
% plot(x(length(M1)+1:end), -M0, 'r')
% plot(x, -[M11 M00], 'b')
% plot(x(length(M11)+1:end), -M00, 'c')
% hold off
% 
% figure(3)
% hold on
% plot(x, -[Q1 Q0], 'g')
% plot(x(length(Q1)+1:end), -Q0, 'r')
% plot(x, -[Q11 Q00], 'b')
% plot(x(length(Q11)+1:end), -Q00, 'c')
% hold off

figure(4)
hold on
plot(x0, Shear_a/(P/t), 'g')
hold off




end