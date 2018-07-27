function [ad1, ad2] = Overlap_Adherent_Load_Distributions(x, t, t_a, F_k, Q_k, M_k, Shr_a, Pl_a, method)
%% Description
% Synt_ax
%   >Adhesive_Stresses_Adherent_Loads(x, F, D_1, A_11, D_0, l_A, l_B, E, t, E_a, G_a, t_a)
% Input
%   >x      : [nxn] matrix with coordinates of the overlap where n is the #
%               of elements in the overlap region
%   >t      : Adherent thickness (symmteric adherents) [m]
%   >t_a    : Adhesive thickness [m]
%   >F_k    : Axial edge load [N/m]
%   >M_k    : Moment edge load [N/m^2]
%   >Q_k    : Shear edge load [N/m]
%   >Shr_a  : Adhesive shear stress [N/m^2]
%   >Pl_a   : Adhesive peel stress [N/m^2]
%   >method : num for numerical or ana for analytical
%
% Output
%   >ad1    : Adherent 1
%               1. Axial load
%               2. Shear load
%               3. Bending moment
%   >ad2    : Adherent 2
%               1. Axial load
%               2. Shear load
%               3. Bending moment
% -------------------------------------------------------------------------
%
%% Code

% Note: Numerical integration instead of the exact integral solutions due
% to the increased computational speed
dx = x(1,end) - x(1,end-1);

switch method
    case 'num'
        % Upper adherent - Axial force
        N_1     = F_k.*triu(ones(size(Shr_a(:,:,1))))-cumsum(Shr_a.*dx,2);
        
        % Upper adherent - Shear force
        Q_1     = Q_k.*triu(ones(size(Pl_a(:,:,1))))-cumsum(Pl_a.*dx,2);
        
        % Upper adherent - Bending moment
        M_1     = M_k.*triu(ones(size(Shr_a(:,:,1))))+cumsum(Q_1*dx,2)-(t+t_a)/2*cumsum(Shr_a.*dx,2);
        
        % Lower adherent - Axial force
        N_2     = cumsum(Shr_a.*dx,2);
        
        % Lower adherent - Shear force
        Q_2     = cumsum(Pl_a.*dx,2);
        
        % Lower adherent - Bending moment
        M_2     = cumsum(Q_2*dx,2)-(t+t_a)/2*cumsum(Shr_a.*dx,2);
    case 'Ana'
        % Upper adherent - Axial force
        Bnd_l   = (F_k*t+6*M_k).*sinh(beta_t*-l_B)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*-l_B;
        Bnd_u   = (F_k*t+6*M_k).*sinh(beta_t*x)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*x;
        N_1     = F_k-(Bnd_u-Bnd_l);
        
        % Upper adherent - Shear force
        Q_1     = Q_k+(-(B_s1+B_s4).*cosh(beta_s.*x).*sin(beta_s.*x)+sinh(beta_s.*x).*cos(beta_s.*x).*(B_s1-B_s4)).*E_a./(t_a.*beta_s)-((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a./(t_a.*beta_s);
        
        % Upper adherent - Bending moment
        M_1     = M_k+Q_k.*x+E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))+(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s)./(t_a.*beta_s)-((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x./(t_a.*beta_s)-Q_k.*x(:,1)-E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))+(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s)./(t_a.*beta_s)+((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x(:,1)./(t_a.*beta_s)-((1./2).*t+(1./2).*t_a).*((3.*F_k.*t-6.*M_k).*x./(8.*t.*l_B)+(F_k.*t+6.*M_k).*sinh(beta_t.*x)./(8.*t.*sinh(beta_t.*l_B))-(3.*F_k.*t-6.*M_k).*x(:,1)./(8.*t.*l_B)-(F_k.*t+6.*M_k).*sinh(beta_t.*x(:,1))./(8.*t.*sinh(beta_t.*l_B)));
        
        % Lower adherent - Axial force
        Bnd_l   = (F_k*t+6*M_k).*sinh(beta_t*-l_B)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*-l_B;
        Bnd_u   = (F_k*t+6*M_k).*sinh(beta_t*x)./(8*t*sinh(beta_t*l_B))+3*(F_k*t-2*M_k)./(8*t*l_B).*x;
        N_2     = Bnd_u-Bnd_l;
        
        % Lower adherent - Shear force
        Q_2     = -(-(B_s1+B_s4).*cosh(beta_s.*x).*sin(beta_s.*x)+sinh(beta_s.*x).*cos(beta_s.*x).*(B_s1-B_s4)).*E_a./(t_a.*beta_s)+((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a./(t_a*beta_s);
        
        % Lower adherent - Bending moment
        M_2     = ((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x./(t_a.*beta_s)-E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))+(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x)+sinh(beta_s.*x)).*cos(beta_s.*x)+(1./4).*sin(beta_s.*x).*(cosh(beta_s.*x)+sinh(beta_s.*x))-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*cos(beta_s.*x)-(1./4).*(cosh(beta_s.*x)-sinh(beta_s.*x)).*sin(beta_s.*x))./beta_s)./(t_a.*beta_s)-((B_s1+B_s4).*cosh(beta_s.*l_B).*sin(beta_s.*l_B)-sinh(beta_s.*l_B).*cos(beta_s.*l_B).*(B_s1-B_s4)).*E_a.*x(:,1)./(t_a.*beta_s)+E_a.*((B_s1-B_s4).*((1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))+(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s+(-B_s1-B_s4).*(-(1./4).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))+(1./4).*sin(beta_s.*x(:,1)).*(cosh(beta_s.*x(:,1))+sinh(beta_s.*x(:,1)))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*cos(beta_s.*x(:,1))-(1./4).*(cosh(beta_s.*x(:,1))-sinh(beta_s.*x(:,1))).*sin(beta_s.*x(:,1)))./beta_s)./(t_a.*beta_s)-((1./2).*t+(1./2).*t_a).*((3.*F_k.*t-6.*M_k).*x./(8.*t.*l_B)+(F_k.*t+6.*M_k).*sinh(beta_t.*x)./(8.*t.*sinh(beta_t.*l_B))-(3.*F_k.*t-6.*M_k).*x(:,1)./(8.*t.*l_B)-(F_k.*t+6.*M_k).*sinh(beta_t.*x(:,1))./(8.*t.*sinh(beta_t.*l_B)));
end

% Store results
ad1.N = N_1;
ad1.Q = Q_1;
ad1.M = M_1;

ad2.N = N_2;
ad2.Q = Q_2;
ad2.M = M_2;
end

