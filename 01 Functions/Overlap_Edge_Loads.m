function [M_k, M_k0, Q_k, Q_k0, V_k, M, Q, w] = Overlap_Edge_Loads(xA, xB, P, D_A, D_B, l_A, l_B, t, t_a, BC)
%% Description
% Syntax
%   >Overlap_Edge_loads(x1, x0, P, D_A, D_B, l_A, l_Btmp, t, ta)
% Input
%   >xA     : [MxN] matrix with coordinates of the free adherent where N
%               is the # of elements in the overlap region and N the # of
%               cracked elements
%   >xB     : [MxN] matrix with coordinates of the overlap where N is the #
%               of elements in the overlap region and M the # of cracked
%               elements
%   >P      : Applied load per unit width along the CLS neutral axis [N/m]
%   >D_A    : Free adherent bending stiffness
%   >D_B    : Overlap bending stiffness (lumped together)
%   >l_A    : Free adherent length [m]
%   >l_B    : Overlap length [m]
%   >t      : Adherent thickness (symmteric adherents) [m]
%   >t_a    : Adhesive thickness [m]
%   >BC     : Support Boundary Conditions
%               1. RR = Roller-Roller (incl. pinned)
%               2. CC = Clamped-Clamped
% Output
%   >M_k    : Moment at x1 = l_A per unit width [N]
%   >M_k0   : Moment at x0 = 0 per unit width [N]
%   >Q_k    : Shear force at x1 = l_A per unit width [N/m]
%   >Q_k0   : Shear force at x0 = 0 per unit width [N/m]
%   >V_k    : Shear force at x1 = l_A per unit width [N/m] on equilibirum 
%               consideration of the overlap region
%   >M      : Structurer with bending moment distribution
%               1) .A for the free adherent
%               2) .B for the overlap region
%   >Q      : Structurer with shear force distribution
%               1) .A for the free adherent
%               2) .B for the overlap region
%   >w      : Structurer with vertical desiplacement distribution
%               1) .A for the free adherent
%               2) .B for the overlap region
%               
% Description
% >Based on Goland and Reissner (1944) solution where both adherents in
%   the overlap region are lumped together
% >Available support boundary conditions:
%       1. Roller-Roller; self derrived following 'The cracked lap shear
%       specimen - A closed form solution' by Lai et al. [1995]
%       2. Clamped-Clamped; self derrived following 'The cracked lap shear
%       specimen - A closed form solution' by Lai et al. [1995]
% -------------------------------------------------------------------------
%
%% Code

% Solution eigenvalues
lambda_A = sqrt(P/D_A);
lambda_B = sqrt(P/D_B);

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t+t_a)./(2*(l_A(1)+l_B(1)));

switch BC
    case 'RR'
        % Reaction forces according to force equilibrium
        R_A = 0;
        M_A = 0;
        
        % Integration constants
        A_1 = 0;
        B_1 = -(sqrt(D_A)*(cosh(lambda_B.*l_B)*t...
            +cosh(lambda_B.*l_B)*t_a+2*alpha*(l_A+l_B)-t-t_a))...
            ./(2*(sqrt(D_A)*cosh(lambda_B.*l_B).*sinh(lambda_A.*l_A)...
            +sqrt(D_B).*cosh(lambda_A.*l_A).*sinh(lambda_B.*l_B)));
        A_0 = (-2*sqrt(D_B/D_A)*B_1.*cosh(lambda_A.*l_A).*sinh(lambda_B.*l_B)...
            -2*alpha*(l_A+l_B)+t+t_a)./(2*cosh(lambda_B.*l_B));
        B_0 = sqrt(D_B/D_A)*B_1.*cosh(lambda_A.*l_A);
        
        % Overlap edge loads and lumped overlap loads
        M_k     = P.*(-A_1.*cosh(lambda_A.*l_A)-B_1.*sinh(lambda_A.*l_A));
        Q_k     = P.*(-A_1.*lambda_A.*sinh(lambda_A.*l_A)-B_1.*lambda_A.*cosh(lambda_A.*l_A));
        M_k0    = P.*(-A_0.*cosh(lambda_B.*0)-B_0.*sinh(lambda_B.*0));
        Q_k0    = P.*(-A_0.*lambda_B.*sinh(lambda_B.*0)-B_0.*lambda_B.*cosh(lambda_B.*0));
        
        % Q_k according to force equilibrium of the undeformed body
        V_k =  (P*(t+t_a)-2*M_k)./(2*l_B);
        
        % Vertical displacement
        w.A = A_1.*cosh(lambda_A.*xA)+B_1.*sinh(lambda_A.*xA)+(alpha+R_A./P).*xA+M_A./P;
        w.B = A_0.*cosh(lambda_B.*xB)+B_0.*sinh(lambda_B.*xB)+alpha.*(l_A+xB)-(t+t_a)/2+(l_A+xB).*R_A./P+M_A./P;
        
        % Free adherent >> CORRECT
        M.A = P.*(-A_1.*cosh(lambda_A.*xA)-B_1.*sinh(lambda_A.*xA));
        Q.A = P.*(-A_1.*lambda_A.*sinh(lambda_A.*xA)-B_1.*lambda_A.*cosh(lambda_A.*xA));
        
        % Overlap region >> CORRECT
        M.B = P.*(-A_0.*cosh(lambda_B.*xB)-B_0.*sinh(lambda_B.*xB));
        Q.B = P.*(-A_0.*lambda_B.*sinh(lambda_B.*xB)-B_0.*lambda_B.*cosh(lambda_B.*xB));
        
    case 'CC'
        %Pre-allocate memory
        R_A = zeros(size(l_B,1),1,size(P,2));
        R_B = zeros(size(l_B,1),1,size(P,2));
        M_A = zeros(size(l_B,1),1,size(P,2));
        M_B = zeros(size(l_B,1),1,size(P,2));
        A_1 = zeros(size(l_B,1),1,size(P,2));
        B_1 = zeros(size(l_B,1),1,size(P,2));
        A_0 = zeros(size(l_B,1),1,size(P,2));
        B_0 = zeros(size(l_B,1),1,size(P,2));
        
        % Solve the 8 equations and 8 unknowns in matrix form
        for i = 1:size(P,3)
            
            % Applied load
            Ptmp        = P(1,1,i);
            lambda_Atmp = lambda_A(1,1,i);
            lambda_Btmp = lambda_B(1,1,i);
            
            for j = 1:size(l_B,1)
                % Free adherent (l_A) and overlap (l_B) length
                l_Btmp = l_B(j);
                l_Atmp = l_A(j);
                
                % Integration constants
                a_mat = [0 0 0.1e1 / Ptmp 0 1 0 0 0; 0.1e1 / Ptmp 0 0 0 0 lambda_Atmp 0 0; 0 0.1e1 / Ptmp * (l_Btmp + l_Atmp) 0.1e1 / Ptmp 0 0 0 cosh(lambda_Btmp * l_Btmp) sinh(lambda_Btmp * l_Btmp); 0 0.1e1 / Ptmp 0 0 0 0 lambda_Btmp * sinh(lambda_Btmp * l_Btmp) lambda_Btmp * cosh(lambda_Btmp * l_Btmp); 0.1e1 / Ptmp * l_Atmp -0.1e1 / Ptmp * l_Atmp 0 0 cosh(lambda_Atmp * l_Atmp) sinh(lambda_Atmp * l_Atmp) -1 0; 0.1e1 / Ptmp -0.1e1 / Ptmp 0 0 lambda_Atmp * sinh(lambda_Atmp * l_Atmp) lambda_Atmp * cosh(lambda_Atmp * l_Atmp) 0 -lambda_Btmp; 1 -1 0 0 0 0 0 0; -l_Btmp - l_Atmp 0 -1 1 0 0 0 0];
                d_mat = [0 -alpha -(l_Btmp + l_Atmp) * alpha + t / 0.2e1 + t_a / 0.2e1 -alpha -t / 0.2e1 - t_a / 0.2e1 0 0 0]';
                b_mat = a_mat\d_mat;
 
                % Isolate integration constants
                R_A(j,1,i) = b_mat(1);
                R_B(j,1,i) = b_mat(2);
                M_A(j,1,i) = b_mat(3);
                M_B(j,1,i) = b_mat(4);
                A_1(j,1,i) = b_mat(5);
                B_1(j,1,i) = b_mat(6);
                A_0(j,1,i) = b_mat(7);
                B_0(j,1,i) = b_mat(8);
            end
        end
        
        % Overlap edge loads and lumped overlap loads
        M_k     = P.*(-A_1.*cosh(lambda_A.*l_A)-B_1.*sinh(lambda_A.*l_A));
        Q_k     = P.*(-A_1.*lambda_A.*sinh(lambda_A.*l_A)-B_1.*lambda_A.*cosh(lambda_A.*l_A));
        M_k0    = P.*(-A_0.*cosh(lambda_B.*0)-B_0.*sinh(lambda_B.*0));
        Q_k0    = P.*(-A_0.*lambda_B.*sinh(lambda_B.*0)-B_0.*lambda_B.*cosh(lambda_B.*0));
        
        % Q_k according to force equilibrium of the undeformed body
        V_k =  (P*(t+t_a)-2*M_k)./(2*l_B);
        
        % Vertical displacement
        w.A = A_1.*cosh(lambda_A.*xA)+B_1.*sinh(lambda_A.*xA)+(alpha+R_A./P).*xA+M_A./P;
        w.B = A_0.*cosh(lambda_B.*xB)+B_0.*sinh(lambda_B.*xB)+alpha.*(l_A+xB)-(t+t_a)/2+(l_A+xB).*R_A./P+M_A./P;
        
        % Free adherent >> CORRECT
        M.A = P.*(-A_1.*cosh(lambda_A.*xA)-B_1.*sinh(lambda_A.*xA));
        Q.A = P.*(-A_1.*lambda_A.*sinh(lambda_A.*xA)-B_1.*lambda_A.*cosh(lambda_A.*xA));
        
        % Overlap region >> CORRECT
        M.B = P.*(-A_0.*cosh(lambda_B.*xB)-B_0.*sinh(lambda_B.*xB));
        Q.B = P.*(-A_0.*lambda_B.*sinh(lambda_B.*xB)-B_0.*lambda_B.*cosh(lambda_B.*xB));
end

end
