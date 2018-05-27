function [Shr_a, Pl_a, Ad_1, Ad_2] = Adhesive_Stresses_LuoTong_vOld2(x, F_k, M_k, Q_k, D_1, A_11, l_B, E, t, E_a, G_a, t_a)
%% Description
% Synt_ax
%   >Overlap_Edge_loads(x, F, D_1, A_11, D_0, l_A, l_B, E, t, E_a, G_a, t_a)
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
%   >Ad_1   : Adherent 1
%               1. Axial load
%               2. Shear load
%               3. Bending moment
%   >Ad_2   : Adherent 2
%               1. Axial load
%               2. Shear load
%               3. Bending moment
% Description
%
% -------------------------------------------------------------------------
%
%% Adhesive Stresses

% From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 33
alpha_k = A_11*t^2/(4*D_1);
alpha_a = (1+alpha_k)/4;

% Import_ant note: applying the k factor for M_k and Q_k from GR
% gives the exact same solution for the adhesive stresses between
% (1) Luo and Tong (2004, 2007) and GR (1944)
%         p = F_k/t;  % From A. Vlot P53
%         u2c = sqrt(3*(1-0.33^2)/2)*l_B/t.*sqrt(p/E); % From A. Vlot P53
%         k = 1./(1+2*sqrt(2)*t_anh(u2c)); % From A. Vlot P53
%         M_k = k.*F_k*t/2; % From A. Vlot P10
%         Q_k = M_k.*sqrt(F_k/D_1); % From A. Vlot P10

% Peel stress integration const_ants (From Luo and Tong 2004)
beta_s = sqrt(2)/2*(24*E_a/(E*t^3*t_a))^(1/4);
beta_t = sqrt(8*G_a/(A_11*t_a));
beta_a = sqrt(alpha_a*beta_t^2);

% Peel stress integration const_ants (From Luo and Tong 2004)
B_s1 =  6/(beta_s^3*E*t^3)*(M_k.*beta_s.*(sinh(beta_s*l_B).*cos(beta_s*l_B)...
    +cosh(beta_s*l_B).*sin(beta_s*l_B))...
    +Q_k.*sinh(beta_s*l_B).*sin(beta_s*l_B))...
    ./(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));
B_s4 =  6/(beta_s^3*E*t^3)*(M_k.*beta_s.*(sinh(beta_s*l_B).*cos(beta_s*l_B)...
    -cosh(beta_s*l_B).*sin(beta_s*l_B))...
    +Q_k.*cosh(beta_s*l_B).*cos(beta_s*l_B))...
    ./(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));

% Adhesive shear stress (Luo and Tong 2007)
Shr_a = beta_t*(F_k*t+6*M_k).*cosh(beta_t*x)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B);
% Adheisve peel stress (Luo and Tong 2004)
Pl_a = 2*E_a/t_a*(B_s1.*sinh(beta_s*x).*sin(beta_s*x)+B_s4.*cosh(beta_s*x).*cos(beta_s*x));

%% Adherent Internal Loads

dx = (x(2)-x(1));

% Upper adherent - Axial force
Bnd_l   = (F_k*t+6*M_k).*sinh(beta_t*-l_B)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*-l_B;
Bnd_u   = (F_k*t+6*M_k).*sinh(beta_t*x)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*x;
N_1     = F_k-(Bnd_u-Bnd_l);
N_1     = F_k-cumsum(Shr_a.*dx,2);

% Upper adherent - Shear force
Q_1     = Q_k+(-(B_s1+B_s4).*cosh(beta_s.*x).*sin(beta_s.*x)+sinh(beta_s.*x).*cos(beta_s.*x).*(B_s1-B_s4)).*E_a./(t_a.*beta_s)-((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a./(t_a.*beta_s);
Q_1     = Q_k-cumsum(Pl_a.*dx,2);

% Upper adherent - Bending moment
M_1  = M_k+cumsum(Q_1*dx,2)-(t+t_a)/2*cumsum(Shr_a.*dx,2);
M_11 = M_k+Q_k.*x+E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))+(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s)./(t_a.*beta_s)-((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x./(t_a.*beta_s)-Q_k.*x(:,1)-E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))+(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s)./(t_a.*beta_s)+((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x(:,1)./(t_a.*beta_s)-((1./2).*t+(1./2).*t_a).*((3.*F_k.*t-6.*M_k).*x./(8.*t.*l_B)+(F_k.*t+6.*M_k).*sinh(beta_t.*x)./(8.*t.*sinh(beta_t.*l_B))-(3.*F_k.*t-6.*M_k).*x(:,1)./(8.*t.*l_B)-(F_k.*t+6.*M_k).*sinh(beta_t.*x(:,1))./(8.*t.*sinh(beta_t.*l_B)));
M_1  = M_1;

% Lower adherent - Axial force
Bnd_l   = (F_k*t+6*M_k).*sinh(beta_t*-l_B)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*-l_B;
Bnd_u   = (F_k*t+6*M_k).*sinh(beta_t*x)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*x;
N_2     = Bnd_u-Bnd_l;
N_2     = cumsum(Shr_a.*dx,2);

% Lower adherent - Shear force
Q_2     = -(-(B_s1+B_s4).*cosh(beta_s.*x).*sin(beta_s.*x)+sinh(beta_s.*x).*cos(beta_s.*x).*(B_s1-B_s4)).*E_a./(t_a.*beta_s)+((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a./(t_a*beta_s);
Q_2     = cumsum(Pl_a.*dx,2);

% Lower adherent - Bending moment
M_2  = cumsum(Q_2*dx,2)-(t+t_a)/2*cumsum(Shr_a.*dx,2);
M_22 = ((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x./(t_a.*beta_s)-E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))+(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s)./(t_a.*beta_s)-((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x(:,1)./(t_a.*beta_s)+E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))+(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s)./(t_a.*beta_s)-((1./2).*t+(1./2).*t_a).*((3.*F_k.*t-6.*M_k).*x./(8.*t.*l_B)+(F_k.*t+6.*M_k).*sinh(beta_t.*x)./(8.*t.*sinh(beta_t.*l_B))-(3.*F_k.*t-6.*M_k).*x(:,1)./(8.*t.*l_B)-(F_k.*t+6.*M_k).*sinh(beta_t.*x(:,1))./(8.*t.*sinh(beta_t.*l_B)));
M_2  = M_22;

% Store results
Ad_1.N = N_1;
Ad_1.Q = Q_1;
Ad_1.M = M_1;

Ad_2.N = N_2;
Ad_2.Q = Q_2;
Ad_2.M = M_2;

end