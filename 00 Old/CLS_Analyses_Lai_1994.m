function [w1, w0, x] = CLS_Analyses_Lai_1994(Load, Strap, Lap, Adhesive, Conditions)
% Author  : Lai et al. (1994)
% Method  : Goland and Reissner
% Remarks : This methods is in error
%           1. Constants a4 and a5 are not the same when taking the
%           compacted or extended equations (appendix A)
%           2. A_1 = correct, B_1 has a insignificant difference w.r.t my
%           own derrivation
%           3. A_0 and B_0 are nowhere near correct.. w.r.t my
%           own derrivation
%           4. It I set t_a = zero in my own derrication, the free adherent
%           has the exact same displacement.
%           5. If I give constants from own derrivation to this solution,
%           results are the same
%           6. Problem is with A_0 and B_0

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

% Define geometry
L   = L_1;
l_a = (L - L_2) + Adhesive.la0;
l_b = L_2 - Adhesive.la0;

D_1 = E_1 * 1/12 * t_1^3;
D_0 = 2 * E_2 * (1/12 * t_1^3 + t_1 * (t_1/2)^2); % >>>> FOR NOW SYMMETRICA

lambda_1 = sqrt(P/D_1);
lambda_0 = sqrt(P/D_0);

% Nondimensional paramaters
P_nd     = sqrt(P / (E_1 * (t_1+t_2)));
q_nd     = l_a/L;
k_nd     = t_1/L;
mu_nd    = t_1/t_2;
S_nd     = E_1/E_2;

% Non-dimensional notations to simplify the unknown constant expressions
%a1 = sinh(q_nd/k_nd * P_nd * sqrt(12 * (1+1/mu_nd)))
a1 = sinh(lambda_1*l_a);

a2 = lambda_1/lambda_0;

%a3 = cosh(q_nd/k_nd * P_nd * sqrt(12 * (1+1/mu_nd)))
a3 = cosh(lambda_1*l_a);

% ERROR HERE
%a4 = sinh((1-q_nd)/k_nd * P_nd * sqrt(12 * (1+1/mu_nd)))
a4 = sinh(lambda_0*l_b);

% ERROR HERE
%a5 = cosh((1-q_nd)/k_nd * P_nd * sqrt(12 * (1+1/mu_nd)))
a5 = cosh(lambda_0*l_b);

a6 = q_nd / (lambda_1*l_a);

a7 = lambda_0 / lambda_1;


% Integration constants
if Conditions.BCs == 1
    % Roller-roller configuration
    delta = (1 + 2 * S_nd * mu_nd + S_nd * mu_nd^2) / (2 * mu_nd * (1 + S_nd * mu_nd));
    
    alpha_nd = (t_1/2 + t_2 - delta * t_1) / L;
    
    R_A_nd = R_A/P;
    M_A_nd = M_A/(P*L);
    
    A_1 = 0;
    B_1 = -alpha_nd*a5 / (a1*a5 + a2*a3*a4);
    A_0 = B_1*a2*a3;
    B_0 = B_1*a1+alpha_nd;
    
elseif Conditions.BCs == 2
elseif Conditions.BCs == 3
end

x1 = 0:l_a/100:l_a;
x0 = 0:l_b/100:l_b;

x = [x1 x0+l_a];

w1 = L * (A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha_nd + R_A_nd) .* x1 / L + M_A_nd);
w0 = L * (A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + (alpha_nd + R_A_nd) .* (x0 + l_a)/L - alpha_nd + M_A_nd);

%figure(1)
%plot(x1,w1,x0+l_a,w0)
end