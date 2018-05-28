function [Minor1, Minor2] = Adherent_Fatigue_Initiation(ABD, PlyStiff, Nx, Mx)
%% Description
% Syntax 
%   >Adherent_Fatigue_Initiation(C, S, N, Q, M)
% Input  
%   >ABD    : Laminate ABD matrix
%   >S      : Laminate ply stiffnesses (Q-matrices)
%   >load   : Load vector [Nx ; Ny ; Nxy ; Mx ; My ; Mxy];
% Output

% Description
%
% -------------------------------------------------------------------------
%
%% Thermal Stresses

%% Mechanical stresses

e   = zeros(size(Nx));
for i = 1:size(Nx,1)
    for j = 1:size(Nx,2)
        
        e = ABD^(-1) * [Nx(i,j,1) ; Ny ; Nxy ; Mx(i,j,1) ; My ; Mxy];
        
        
        
    end
end

Minor1 = 1;
Minor2 = 1;