function [Shr_a, Pl_a] = Overlap_Adhesive_Stresses(x, F_k, M_k, Q_k, l_B, E, t, E_a, G_a, t_a)
%% Description
% Synt_ax
%   >Overlap_Adhesive_Stresses(x, F_k, M_k, Q_k, l_B, E, t, E_a, G_a, t_a)
% Input
%   >x      : [nxn] matrix with coordinates of the overlap where n is the #
%               of elements in the overlap region
%   >F_k    : Axial edge load [N/m]
%   >M_k    : Moment edge load [N/m^2]
%   >Q_k    : Shear edge load [N/m]
%   >D_1    : Free adherent bending stiffness
%   >A_11   : Free adherent axial stiffness
%   >D_0    : Overlap bending stiffness (lumped together)
%   >l_B    : Overlap length [m]
%   >E      : Adherent Young's Modulus (symmteric adherents) [Pa]
%   >t      : Adherent thickness (symmteric adherents) [m]
%   >E_a    : Adhesive Young's Modulus [Pa]
%   >G_a    : Adhesive Shear Modulus [Pa]
%   >t_a    : Adhesive thickness [m]
%
% Output
%   >Shr_a  : Adhesive shear stress [N/m^2]
%   >Pl_a   : Adhesive peel stress [N/m^2]
%
% Description
%
% -------------------------------------------------------------------------
%
% Important note: applying the k factor for M_k and Q_k from GR (from Vlot)
% gives the exact same solution for the adhesive stresses between
% (1) Luo and Tong (2004, 2007) and GR (1944)
% p = F_k/t;  % From A. Vlot P53
% u2c = sqrt(3*(1-0.33^2)/2)*l_B/t.*sqrt(p/E);   % From A. Vlot P53
% k = 1./(1+2*sqrt(2)*tanh(u2c));                % From A. Vlot P53
% M_k = k.*F_k*t/2;                              % From A. Vlot P10
% Q_k = -M_k.*sqrt(F_k/D_1);                     % From A. Vlot P10 (- sign instead of + sign)

% Peel stress integration constants 
%   > Source: Modeling of Adhesively Bonded Joints, page 33
beta_s = sqrt(2)/2*(24*E_a/(E*t^3*t_a))^(1/4);
beta_t = sqrt(8*G_a/(E*t*t_a));

% Integration constants
%   > Source: Modeling of Adhesively Bonded Joints, page 42
B_s1 = M_k.*(sinh(beta_s*l_B).*cos(beta_s*l_B)+cosh(beta_s*l_B).*sin(beta_s*l_B))...
    +Q_k/beta_s.*sinh(beta_s*l_B).*sin(beta_s*l_B);
B_s4 = M_k.*(sinh(beta_s*l_B).*cos(beta_s*l_B)-cosh(beta_s*l_B).*sin(beta_s*l_B))...
    +Q_k/beta_s.*cosh(beta_s*l_B).*cos(beta_s*l_B);

% Adhesive shear stress
%   > Source: Modeling of Adhesively Bonded Joints, page 41
Shr_a = beta_t*(F_k*t+6*M_k).*cosh(beta_t*x)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B);
%Shr_a = Shr_a.*triu(ones(size(Shr_a(:,:,1))));

% Adhesive peel stress
%   > Source: Modeling of Adhesively Bonded Joints, page 41
Pl_a = 2*beta_s^2*(B_s1.*sinh(beta_s*x).*sin(beta_s*x)+B_s4.*cosh(beta_s*x).*cos(beta_s*x))...
    ./(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));
%Pl_a = Pl_a.*triu(ones(size(Pl_a(:,:,1))));

end