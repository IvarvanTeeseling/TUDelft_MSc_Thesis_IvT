function [Shr_a, Pl_a] = Overlap_Adhesive_Stresses_LT(x, F_k, M_k, Q_k, l_B, E, G, t, E_a, G_a, t_a, EA, EI)
%% Description
% Synt_ax
%   >Overlap_Adhesive_Stresses_LT(x, F_k, M_k, Q_k, l_B, E, t, E_a, G_a, t_a)
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
%% Code

% Shear Stress
k_t = G_a/t_a;

alpha_k1 = EA*t^2/(4*EI);
alpha_k2 = EA*t^2/(4*EI);

alpha_a1 = (1+alpha_k1)/4;
alpha_a2 = (1+alpha_k2)/4;

k_1d = 4*alpha_a1/EA+4*alpha_a2/EA;

beta_tt = sqrt(k_1d);

beta_t = sqrt(k_t)*beta_tt;

H_nd = -F_k/EA-1/2*t*M_k/EI;

H_n = k_t*H_nd;

Shr_a = -H_n/beta_t.*exp(-beta_t*x);

% Peel stress
k_4d = 1/EI+1/EI;

k_5d = 1/(G*t)+1/(G*t);

C_f = 4*k_4d/(E_a*k_5d^2);

if t_a>1/C_f
    k_s = E_a/t_a;
    
    beta_tg = sqrt(k_5d/4);
    
    beta_t1 = beta_tg*sqrt(sqrt(C_f*t_a)+1);
    beta_s1 = sqrt(k_s)*beta_t1;
    
    beta_t2 = beta_tg*sqrt(sqrt(C_f*t_a)-1);
    beta_s2 = sqrt(k_s)*beta_t2;
    
    H_md = M_k/EI;
    
    H_qd = Q_k/EI;

    beta_ts = (k_4d/(4*k_s))^(1/4);
    
    B_3 = -(sqrt(sqrt(C_f*t_a)+1)*H_md+(beta_tg*H_qd)/(beta_ts^2*sqrt(k_s)))...
        /(2*beta_ts^2*sqrt(sqrt(C_f*t_a)-1));
    
    B_4 = 1/(2*beta_ts^2)*(H_md+(beta_tg*sqrt(sqrt(C_f*t_a)+1))/(beta_ts^2*sqrt(k_s))*H_qd);
    
    Pl_a = exp(-beta_s1*x).*(B_3.*sin(beta_s2*x)+B_4.*cos(beta_s2*x));
    
elseif t_a<C_f
    
end

end
