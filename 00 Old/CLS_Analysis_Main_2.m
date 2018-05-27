clear all; close all; clc

%%
% >>>> Input
% Applied load [N/m]
P   = 20000; 
% Applied load Y direction in A [N/m]
R_A = 0;
% Applied moment in A [N]
M_A = 0;

% Strap adherent modulus [Pa]
E_1  = 72.45e9;     
% Strap adherent thicness [m]
t_1  = 0.00318;
% Strap adherent poisson [-]
v_1  = 0.33; 
% Strap length [m]
L_1  = 0.405;

% Lap adherent modulus [Pa]
E_2  = 72.45e9; 
% Lap adherent thicness [m]
t_2  = 0.00318;        
% Lap adherent poisson [-]
v_2  = 0.33;         
% Lap length [m]
L_2  = 0.254;

% Adhesive modulus [Pa]
E_a  = 1.932e9; 
% Adhesive  thicness [m]
t_a  = 0.00013;        
% Adhesive poisson [-]
v_a  = 0.4;  

% Initial crack length [m]
la0  = 0;

% Stress conditions
PlaneStrain = 0;    % 1 = yes, 0 = no
BCs         = 1;    % 1 = r-r, 2 = r-c, 3 = c-c.
                    % Note: left is free adherent, right is overlap

% >>>> Correct modulus to plane strain if needed
if PlaneStrain == 1
   E_1 = E_1 / (1-v_1^2);
   E_2 = E_2 / (1-v_2^2);
end  

% >>>> Part 1: Forces and vertical displacement functions

l_a  = (L_1 - L_2) + la0;
l_b  = L_2 - la0;
L   = L_1;

%% 
% Author  : Lai et al. (1994)
% Method  : Goland and Reissner
% Remarks : Stresses added as they are not in the paper

D_1 = E_1 * 1/12 * t_1^3;
D_0 = 2 * (D_1 + t_1 * (t_1/2)^2); % >>>> FOR NOW


%% Modeling of Adhesively bonded Joints

c = l_b;
l = l_a;

alpha = (t_2 + t_a)/(2*(l+c));

beta_k = sqrt(P/D_1);
beta_0 = sqrt(P/D_0);

A1 = -1 * (beta_0*(t_2 + t_a) * cosh(beta_0*c)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
B1 = -1 * (beta_k*(t_2 + t_a) * cosh(beta_k*l)) / (2*(beta_0*sinh(beta_k*l)*cosh(beta_0*c) + beta_k*cosh(beta_k*l)*sinh(beta_0*c)));
A0 = 0;
B0 = 0;

x1 = 0:l_a/100:l_a;
x0 = -l_b:l_b/100:0;

w1 = A1 * sinh(beta_k * x1) + alpha * x1;
w0 = B1 * sinh(beta_0 * x0) + alpha * x0;

plot(x1,w1,'--',x0+l_a+l_b,w0,':')