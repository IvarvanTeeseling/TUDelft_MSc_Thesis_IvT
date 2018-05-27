function [serr, mr] = Calc_StrainEnergyReleaseRate(P, M_k, M_0, E, t, EA, D, method)
%% Description
% Syntax
%   >[serr, mr] = Calc_StrainEnergyReleaseRate(p, mk, E, t, AE, D)
% Input
%   >F      : Tensile load (left of the crack tip)
%   >M_k    : Overlap edge bending moment (left of the crack tip)
%   >M_0    : Overlap bending moment (right of the crack tip)
%   >E      : Adherent Young's Modulus
%   >P      : Applied load per unit width along the CLS neutral axis [N/m]
%   >D      : Free adherent bending stiffness
%   >t      : Adherent thickness (symmteric adherents) [m]
%   >method : SERR calculation method
%               1. Verreman - On the fracture paramaters in a clamped
%               cracked lap shear adhesive joint [1992]
%               2. Lai - The cracked lap shear specimen revisited -
%               A closed form solution [1996]
%               3. Fernlund1991 - Failure load prediction of structural
%               adhesive joints - Part 1 [1991]
%               4. Brussat - Fracture mechanics for structural adhesive
%               bonds - Final report [1977]
%               5. Fernlund1994 - Fracture load predictions for adhesive
%               joints [1994]
% Output
%   >serr   : Strain Energy Release Rate
%               1. serr.G for the total SERR
%               2. serr.GI for the mode I SERR
%               3. serr.GII for the mode II SERR
%   >mr     : SERR Mode Ratio (GII/G)
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
switch method
    case 'Verreman'
        % Adhesive stresses to get Mode I and Mode II
        GI = ta/(2*Ea)*peel_a.^2;
        GII = ta/(2*Ga)*shear_a.^2;
        G = GI+GII;
        mr = GII./G;
    case 'Lai'
        % Used Sue and Hutchinson 'Interface crack between two elastic
        % layers' [1990] to find (1) total SERR and (2) the mode
        % partitioning using the Suo and Hutchinson (1990) approach
        mu = 1;
        sm = 1;
        delta = 1;
        I = sm*((delta-1/mu)^2-(delta-1/mu)+1/3)+delta/mu*(delta-1/mu)+1/(3*mu^3);
        A = 1/mu + sm;
        
        G = 1/(2*E)*((P.^2/t+12*M_k.^2./t^3)+(-P.^2/(A*t)-M_0.^2/(I*t^3)));
        GI = 0.75*G;
        GII = 0.25*G;
        mr = GII./G;
        % TO DO: include mode ratio from Suo & Hutchinson
    case 'Fern1und1991'
        % Used the J-intregral and characteristics of similar adherends to
        % find Mode I and Mode II
        GI = M_k.^2/(4*D);
        GII = P.^2/(4*EA)+(4*M_k.^2-M_0.^2)./(16*D);
        G = GI+GII;
        mr = GII./G;
    case 'Brussat'
        % Infintely long specimen where SERR is found using energy balance
        % considerations and the MR is an approximation
        G = P.^2/(4*E*t);
        GI = 0.75*G;
        GII = 0.25*G;
        mr = GII./G;
    case 'Fernlund1994'
        % Used the J-integral for the total SERR and Sue and Hutchinson 'Interface crack between two elastic
        % layers' [1990] for the mode partitioning
end

% Set output
serr.G = G;
serr.GI = GI;
serr.GII = GII;

end