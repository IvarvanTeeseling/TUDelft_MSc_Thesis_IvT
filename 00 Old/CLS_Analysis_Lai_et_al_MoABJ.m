function [w1, w0, M1, M0, x, shear_a, peel_a] = CLS_Analysis_MoABJ(Load, Strap, Lap, Adhesive, Conditions)
% Paper     > Modeling of Adhesively Bonded Joints
% Author    > Lucas Filipe Martins da Silva (2008)
% Method    > Based on Lai et al (1994)
% Includes  > Normalized solutions for (1) overlap loads and
%           displacements, (2) SERR, (3) Mode Ratio and (4) Fracture
%           Efficieny
%           > Solutions for the BCs (1) C-C, (2) R-C, (3) R-R
%
% Remarks   > Adhesive stresses added as they are not included
%           > Adherent stresses added as they are not included
%           > Integration constants self derived; errors found in their
%           paper
%           > Compared and confirmed with solution from the self derrived
%           solution of Lai et al. (1994).
%% Load Input Data

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

% Re-Define to match the definitions used in the book
l = l_a;
c = l_b;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t_2 + t_a)/(2*(l+c));

% Bending stiffness - free adherent
D_11 = E_1 * 1/12 * t_1^3;
% Bending stiffness - overlap region (lumped together)
% To Do: include assymetric cross-sections
% To Do: include t_a
D_00 = 2 * E_2 * (1/12 * t_1^3 + t_1 * (t_1/2)^2);

% Axial stiffness - free adherent
A_11 = E_1 * t_1;
% Axial stiffness - overlap region (lumped together)
% Note: adhesive ignored
A_22 = E_1 * (t_1 + t_2);

%% Applied load

% Note: F in this solution is the horizontal component of the load P (P =
% applied along the neutral line axis)
F = P * (l+c)/ sqrt((l+c)^2+(t_a/2+t_2/2)^2);

% To compare with plots in the book
%F = (8/c)^2*D_11;

%% Axis preperation

% x-vector spanning the free adherent length 
x1 = 0:l_a/100:l_a;
% x-vector spanning the overlap length
x0 = -l_b:l_b/100:0;
% x-vector spanning the entire length
x = [x1 x0+l_b+l_a];

%% Step 1 - Displacements and forces

% Roller-Roller BC
beta_k = sqrt(F/D_11);
beta_0 = sqrt(F/D_00);

% Y instead of A for the integration constants
Y_1 = -1 * (beta_0*(t_2 + t_a) * cosh(beta_0*c)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
B_1 = -1 * (beta_k*(t_2 + t_a) * cosh(beta_k*l)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
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

% Clamped-Clamped BC
%A_1 = beta_0*(2*cosh(beta_0*c)*sinh(beta_0*c)*beta_0*t_2+2*cosh(beta_0*c)*sinh(beta_0*c)*beta_0*t_a+cosh(beta_0*c)*sinh(beta_k*l)*beta_k*t_2+cosh(beta_0*c)*sinh(beta_k*l)*beta_k*t_a-cosh(beta_k*l)*sinh(beta_0*c)*beta_0*t_2-cosh(beta_k*l)*sinh(beta_0*c)*beta_0*t_a-beta_0*sinh(beta_0*c)*t_2-beta_0*sinh(beta_0*c)*t_a-sinh(beta_k*l)*beta_k*t_2-sinh(beta_k*l)*beta_k*t_a)/(2*(cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*c+cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*l+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*c+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*l-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*c+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_0*c)^2*beta_0*beta_k-2*cosh(beta_0*c)*cosh(beta_k*l)*beta_0*beta_k+cosh(beta_k*l)^2*beta_0*beta_k-sinh(beta_0*c)^2*beta_0*beta_k-sinh(beta_0*c)*sinh(beta_k*l)*beta_0^2-sinh(beta_0*c)*sinh(beta_k*l)*beta_k^2-sinh(beta_k*l)^2*beta_0*beta_k));
%A_2 = (-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c*t_2-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c*t_a-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l*t_2-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l*t_a-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*c*t_2-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*c*t_a-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*l*t_2-cosh(beta_0*c)*sinh(beta_0*c)*beta_k*l*t_a+cosh(beta_0*c)^2*beta_k*t_2+cosh(beta_0*c)^2*beta_k*t_a-cosh(beta_0*c)*cosh(beta_k*l)*beta_k*t_2-cosh(beta_0*c)*cosh(beta_k*l)*beta_k*t_a+sinh(beta_0*c)^2*beta_k*t_2+sinh(beta_0*c)^2*beta_k*t_a+sinh(beta_0*c)*sinh(beta_k*l)*beta_0*t_2+sinh(beta_0*c)*sinh(beta_k*l)*beta_0*t_a+sinh(beta_0*c)*beta_k*c*t_2+sinh(beta_0*c)*beta_k*c*t_a+sinh(beta_0*c)*beta_k*l*t_2+sinh(beta_0*c)*beta_k*l*t_a-cosh(beta_0*c)*beta_k*t_2-cosh(beta_0*c)*beta_k*t_a+cosh(beta_k*l)*beta_k*t_2+cosh(beta_k*l)*beta_k*t_a)*beta_0/(2*(cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*c+cosh(beta_0*c)*sinh(beta_0*c)*beta_0^2*beta_k*l+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*c+cosh(beta_0*c)*sinh(beta_k*l)*beta_0*beta_k^2*l-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*c-cosh(beta_0*c)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*c+cosh(beta_k*l)*sinh(beta_0*c)*beta_0*beta_k*l+cosh(beta_0*c)^2*beta_0*beta_k-2*cosh(beta_0*c)*cosh(beta_k*l)*beta_0*beta_k+cosh(beta_k*l)^2*beta_0*beta_k-sinh(beta_0*c)^2*beta_0*beta_k-sinh(beta_0*c)*sinh(beta_k*l)*beta_0^2-sinh(beta_0*c)*sinh(beta_k*l)*beta_k^2-sinh(beta_k*l)^2*beta_0*beta_k));

Avoids the MATLAB rounding error
a11 = beta_k*beta_0*(l+c)*cosh(beta_0*c)-beta_k*sinh(beta_0*c)-beta_0*sinh(beta_k*l);
a12 = beta_0*(cosh(beta_0*c)-cosh(beta_k*l));
c1  = (t_2+t_a)/2*beta_0*(cosh(beta_0*c)-1);

a21 = beta_k*(cosh(beta_0*c)-cosh(beta_k*l))-beta_k*(l+c)*sinh(beta_0*c);
a22 = -1*(beta_k*sinh(beta_k*l)+beta_0*sinh(beta_0*c));
c2  = (t_2+t_a)/2*beta_0*sinh(beta_0*c);

Values =  [a11 a12 ; a21 a22]^(-1) * [c1 ; c2];

A_1 = Values(1);
A_2 = Values(2);

B_1 = beta_k/beta_0*A_1;
B_2 = A_2+A_1*beta_k*(l+c)+(t_2+t_a)/2;

M_A = -A_2*F;
R_A = -beta_k*A_1*F;
M_B = -B_2*F;
R_B = beta_0*B_1*F;

M_B - M_A + F*(t_2 + t_a)/2 - R_A*(l+c)

w1 = A_1 * sinh(beta_k*x1) + A_2*cosh(beta_k*x1) + R_A/F*x1 + M_A/F;
w0 = B_1 * sinh(beta_0*x0) + B_2*cosh(beta_0*x0) - R_B/F*x0 + M_B/F;

%w11 = A_1*sqrt(F)*sinh(sqrt(F)*x_1/sqrt(D_11))/sqrt(D_11)+B_1*sqrt(F)*cosh(sqrt(F)*x_1/sqrt(D_11))/sqrt(D_11)-R_A/F;
%w22 = A_0*sqrt(F)*sinh(sqrt(F)*x_0/sqrt(D_2))/sqrt(D_2)+B_0*sqrt(F)*cosh(sqrt(F)*x_0/sqrt(D_2))/sqrt(D_2)+R_A/F;

figure(1)
hold on
plot(x, -[w1 w0]/t_1, 'g')
plot(x(length(w1)+1:end), -w0/t_1, 'r')
hold off

% w11 = A_1*(sinh(beta_k*x1)-beta_k*x1)+A_2*(cosh(beta_k*x1)-1);
% w00 = A_1*(beta_k/beta_0*(sinh(beta_0*x0)-beta_0*x0)+ ...
%             beta_k*(l+c)*(cosh(beta_0*x0)-1)) ...
%             +A_2*(cosh(beta_0*x0)-1)+(t_2+t_a)/2*(cosh(beta_0*x0)-1);
       
        
% Edge Moment Factor
%k = 1/(1+(beta_k/beta_0)*tanh(beta_0*c)*coth(beta_k*l));
%k = 1/(1+sqrt(D_00/D_11)*tanh(sqrt(D_11/D_00)*beta_k*c)*coth(beta_k*l));

%% Step 2 - Stresses

% Find the integration constants
alpha_k = A_11 * t_1^2 / (4 * D_11);
alpha_a = (1 + alpha_k) / 4;
 
M_k = M1(end);
 
V_k = 0;
 
beta_s = sqrt(2)/2 * (2*E_a / (D_11*t_a))^(1/4);
beta_t = sqrt(8 * G_a / (A_11*t_a));
beta_a = sqrt(alpha_a * beta_t^2);

B_s1 = M_k*(sinh(beta_s*c*cos(beta_s*c))+cosh(beta_s*c)*sin(beta_0*c)) + V_k/beta_s*sinh(beta_s*c)*sin(beta_s*c);
B_s4 = M_k*(sinh(beta_s*c*cos(beta_s*c))-cosh(beta_s*c)*sin(beta_0*c)) + V_k/beta_s*cosh(beta_s*c)*cos(beta_s*c);

% Adhesive shear stress
shear_a = beta_t*(F*t_1+6*M_k)*cosh(beta_t*x0)/(8*t_1*sinh(beta_t*c) + 3*(F*t_1-2*M_k)/(8*t_1*c));
% Adheisve peel stress
peel_a = 2*beta_s^2*(B_s1*sinh(beta_s*x0).*sin(beta_s*x0)+B_s4*cosh(beta_s*x0).*cos(beta_s*x0)) / (sinh(2*beta_s*c)+sin(2*beta_s));


end