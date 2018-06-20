function [M, Q, V] = Overlap_Edge_Loads(x1, x0, P, D_1, D_0, l_A, l_B, t, t_a, BC)
%% Description
% Syntax
%   >Overlap_Edge_loads(x1, x0, P, D_1, D_0, l_A, l_Btmp, t, ta)
% Input
%   >x1     : [nxn] matrix with coordinates of the free adherent where n
%               is the # of elements in the overlap region
%   >x0     : [nxn] matrix with coordinates of the overlap where n is the #
%               of elements in the overlap region
%   >P      : Applied load per unit width along the CLS neutral axis [N/m]
%   >D_1    : Free adherent bending stiffness
%   >D_0    : Overlap bending stiffness (lumped together)
%   >l_A    : Free adherent length [m]
%   >l_B    : Overlap length [m]
%   >t      : Adherent thickness (symmteric adherents) [m]
%   >t_a    : Adhesive thickness [m]
%   >BC     : Support Boundary Conditions
%               1. RR = Roller-Roller (incl. pinned)
%               2. CC = Clamped-Clamped
% Output
%   >M      : Overlap edge moment per unit width [N/m^2]
%               1. M.A for the moment at x1 = l_A
%               2. M.B for the moment at x0 = 0
%   >Q      : Overlap edge load per unit width w.r.t deformed axis [N/m]
%               1. Q.A for the shear force at x1 = l_A
%               2. Q.B for the shear force at x0 = 0
%   >V      : Overlap edge shear load based on equilibirum consideration 
%               of the overlap region
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
lambda_1 = sqrt(P/D_1);
lambda_0 = sqrt(P/D_0);

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t + t_a)./(2*(l_A(1)+l_B(1)));

switch BC
    case 'RR'
        
        % Reaction forces according to force equilibrium
        R_A = 0;
        M_A = 0;
        
        % Integration constants
        A_1 =   0;
        B_1 =   -(sqrt(D_1)*(cosh(lambda_0.*l_B)*t...
            +cosh(lambda_0.*l_B)*t_a+2*alpha*(l_A+l_B)-t-t_a))...
            ./(2*(sqrt(D_1)*cosh(lambda_0.*l_B).*sinh(lambda_1.*l_A)...
            +sqrt(D_0).*cosh(lambda_1.*l_A).*sinh(lambda_0.*l_B)));
        A_0 =   (-2*sqrt(D_0/D_1)*B_1.*cosh(lambda_1.*l_A).*sinh(lambda_0.*l_B)...
            -2*alpha*(l_A+l_B)+t+t_a)./(2*cosh(lambda_0.*l_B));
        B_0 =   sqrt(D_0/D_1)*B_1.*cosh(lambda_1.*l_A);
        
        % Overlap edge loads
        M_k = P.*(-A_1.*cosh(lambda_1.*l_A)-B_1.*sinh(lambda_1.*l_A));
        Q_k = P.*(-A_1.*lambda_1.*sinh(lambda_1.*l_A)-B_1.*lambda_1.*sinh(lambda_1.*l_A));
        M_0 = P.*(-A_0.*cosh(lambda_0.*0)-B_0.*sinh(lambda_0.*0));
        Q_0 = P.*(-A_0.*lambda_0.*sinh(lambda_0.*0)-B_0.*lambda_0.*sinh(lambda_0.*0));
        
        % Q_k according to force equilibrium of the undeformed body
        V_k =  (P*(t+t_a)-2*M_k)./(2*l_B);
        
        % Set output
        M.A = M_k;
        M.B = M_0;
        Q.A = Q_k;
        Q.B = Q_0;
        V   = V_k;
        
        % Vertical displacement
        %w1 = A_1.*cosh(lambda_1.*x1)+B_1.*sinh(lambda_1.*x1)+(alpha+R_A./P).*x1+M_A./P;
        %w0 = A_0.*cosh(lambda_0.*x0)+B_0.*sinh(lambda_0.*x0)+alpha.*(l_A+x0)-(t+ta)/2+(l_A+x0).*R_A./P+M_A./P;
        
        % Bending moment in the free adherent >> CORRECT
        %M1 = P.*(-A_1.*cosh(lambda_1.*x1)-B_1.*sinh(lambda_1.*x1));
        %Q1 = P.*(-A_1.*lambda_1.*cosh(lambda_1.*x1)-B_1.*lambda_1.*sinh(lambda_1.*x1));
        
        % Bending moment in the overlap region >> CORRECT
        %M0 = P.*(-A_0.*cosh(lambda_0.*x0)-B_0.*sinh(lambda_0.*x0));
        %Q0 = P.*(-A_0.*lambda_0.*sinh(lambda_0.*x0)-B_0.*lambda_0.*sinh(lambda_0.*x0));   
        
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
            lambda_1tmp = lambda_1(1,1,i);
            lambda_0tmp = lambda_0(1,1,i);
            
            for j = 1:size(l_B,1)
                % Free adherent (l_A) and overlap (l_B) length
                l_Btmp = l_B(j);
                l_Atmp = l_A(j);
                
                % Integration constants
                a_mat = [0 0 0.1e1 / Ptmp 0 1 0 0 0; 0.1e1 / Ptmp 0 0 0 0 lambda_1tmp 0 0; 0 0.1e1 / Ptmp * (l_Btmp + l_Atmp) 0.1e1 / Ptmp 0 0 0 cosh(lambda_0tmp * l_Btmp) sinh(lambda_0tmp * l_Btmp); 0 0.1e1 / Ptmp 0 0 0 0 lambda_0tmp * sinh(lambda_0tmp * l_Btmp) lambda_0tmp * cosh(lambda_0tmp * l_Btmp); 0.1e1 / Ptmp * l_Atmp -0.1e1 / Ptmp * l_Atmp 0 0 cosh(lambda_1tmp * l_Atmp) sinh(lambda_1tmp * l_Atmp) -1 0; 0.1e1 / Ptmp -0.1e1 / Ptmp 0 0 lambda_1tmp * sinh(lambda_1tmp * l_Atmp) lambda_1tmp * cosh(lambda_1tmp * l_Atmp) 0 -lambda_0tmp; 1 -1 0 0 0 0 0 0; -l_Btmp - l_Atmp 0 -1 1 0 0 0 0];
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
        M_k = P.*(-A_1.*cosh(lambda_1.*l_A)-B_1.*sinh(lambda_1.*l_A));
        Q_k = P.*(-A_1.*lambda_1.*sinh(lambda_1.*l_A)-B_1.*lambda_1.*sinh(lambda_1.*l_A));
        M_0 = P.*(-A_0.*cosh(lambda_0.*0)-B_0.*sinh(lambda_0.*0));
        Q_0 = P.*(-A_0.*lambda_0.*sinh(lambda_0.*0)-B_0.*lambda_0.*sinh(lambda_0.*0));
        
        % Q_k according to force equilibrium of the undeformed body
        V_k =  (P*(t+t_a)-2*M_k)./(2*l_B);
        
        % Set output
        M.A = M_k;
        M.B = M_0;
        Q.A = Q_k;
        Q.B = Q_0;
        V   = V_k;
        
        % Vertical displacement
        %w1 = A_1.*cosh(lambda_1.*x1)+B_1.*sinh(lambda_1.*x1)+(alpha+R_A./P).*x1+M_A./P;
        %w0 = A_0.*cosh(lambda_0.*x0)+B_0.*sinh(lambda_0.*x0)+alpha.*(l_A+x0)-(t+ta)/2+(l_A+x0).*R_A./P+M_A./P;
        
        % Bending moment in the free adherent >> CORRECT
        %M1 = P.*(-A_1.*cosh(lambda_1.*x1)-B_1.*sinh(lambda_1.*x1));
        %Q1 = P.*(-A_1.*lambda_1.*cosh(lambda_1.*x1)-B_1.*lambda_1.*sinh(lambda_1.*x1));
        
        % Bending moment in the overlap region >> CORRECT
        %M0 = P.*(-A_0.*cosh(lambda_0.*x0)-B_0.*sinh(lambda_0.*x0));
        %Q0 = P.*(-A_0.*lambda_0.*sinh(lambda_0.*x0)-B_0.*lambda_0.*sinh(lambda_0.*x0));      
end

end
