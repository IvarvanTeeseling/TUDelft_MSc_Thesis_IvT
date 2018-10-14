function [serr] = Strain_Energy_Release_Rate(varargin)
%% Description
% Syntax
%   >[serr, mr] = Strain_Energy_Release_Rate(varargin)
% Input
%   >varargin (see definition below per method)
% Output
%   >serr   : Strain Energy Release Rate
%               1. serr.G for the total SERR
%               2. serr.GI for the mode I SERR
%               3. serr.GII for the mode II SERR
%
% Description
%
% -------------------------------------------------------------------------
%
%% Input Check

mstr = {'Verreman1992' 'Lai1996' 'Fern1und1991' 'Fern1und1994' 'Brussat1977'};

if any(strcmp(varargin{1},mstr))
    method = varargin{1};
else
    error('Error. Please enter a valid method.')
end

switch method
    case 'Verreman1992'
        if length(varargin) == 6
            t_a     = varargin{2};
            E_a     = varargin{3};
            G_a     = varargin{4};
            Peel_a  = varargin{5};
            Shear_a = varargin{6};
        else
            error(['Error. Method <' method '> requires 6 inputs: method, t_a, E_a, G_a, Peel_a, Shear_a. You gave me ' num2str(length(varargin)) ' inputs.'])
        end
    case 'Lai1996'
        if length(varargin) == 9
            E       = varargin{2};
            t       = varargin{3};
            EA_1    = varargin{4};
            EI_1    = varargin{5};
            EI_0    = varargin{6};
            P       = varargin{7};
            M_k     = varargin{8};
            M_0     = varargin{9};
        else
            error(['Error. Method <' method '> requires 9 inputs: method, E, t, EA_1, EI_1, EI_0, P, M_k, M_0. You gave me ' num2str(length(varargin)) ' inputs.'])
        end
    case 'Fern1und1991'
        if length(varargin) == 6
            EA_1    = varargin{2};
            EI_1    = varargin{3};
            P       = varargin{4};
            M_k     = varargin{5};
            M_0     = varargin{6};
        else
            error(['Error. Method <' method '> requires 6 inputs: method, EA_1, EI_1, P, M_k, M_0. You gave me ' num2str(length(varargin)) ' inputs.'])
        end
    case 'Fern1und1994'
        if length(varargin) == 8
            t       = varargin{2};
            EA_1    = varargin{3};
            EI_1    = varargin{4};
            EI_0    = varargin{5};
            P       = varargin{6};
            M_k     = varargin{7};
            M_0     = varargin{8};
        else
            error(['Error. Method <' method '> requires 8 inputs: method, t, EA_1, EI_1, EI_0, P, M_k, M_0. You gave me ' num2str(length(varargin)) ' inputs.'])
        end
    case 'Brussat1977'
        if length(varargin) == 9
            t       = varargin{2};
            E       = varargin{3};
            EA_1    = varargin{4};
            EI_1    = varargin{5};
            EI_0    = varargin{6};
            P       = varargin{7};
            M_k     = varargin{8};
            M_0     = varargin{9};
        else
            error(['Error. Method <' method '> requires 9 inputs: method, t, E, EA_1, EI_1, EI_0, P, M_k, M_0. You gave me ' num2str(length(varargin)) ' inputs.'])
        end
end

%% Code

switch method
    case 'Verreman1992'
        % Source:
        %   > "On the fracture parameters in a clamped cracked lap shear 
        %   adhesive joint" (1992) by Edde, F. and Verreman, Y. 
        %
        % Description:
        %   > SERR Mode I and Mode II using the adhesive stress formulation
        
        GI  = t_a/(2*E_a)*Peel_a.^2;
        GII = t_a/(2*G_a)*Shear_a.^2;
        G   = GI+GII;
        mr  = GII./G;
    case 'Lai1996'
        %   > NOT FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % Source:
        %   > "The cracked lap shear specimen revisited - A closed 
        %   form solution" (1996) by Lai et al.
        %
        % Description:
        %   > Suo and Hutchinson to find (1) the total SERR and (2) perform
        %   the mode partitioning
        
        mu      = 1;
        sm      = 1;
        delta   = 1;
        I       = sm*((delta-1/mu)^2-(delta-1/mu)+1/3)+delta/mu*(delta-1/mu)+1/(3*mu^3);
        A       = 1/mu + sm;
        
        G   = 1/(2*E)*((P.^2/t+12*M_k.^2./t^3)+(-P.^2/(A*t)-M_0.^2/(I*t^3)));
        GI  = 0.75*G;
        GII = 0.25*G;

    case 'Fern1und1991'
        % Source:
        %   > "Failure load prediction of structural adhesive joints -
        %   Part 1" (1991) by Fernlund, G. Spelt, J.K.
        %
        % Description:
        %   > SERR for the J-integral and mode partitioning only valid for
        %   symmetric, balanced CLS joints
        
        GI      = M_k.^2/(4*EI_1);
        G_II_f  = P.^2/(4*EA_1);
        G_II_b  = (4*M_k.^2-M_0.^2)./(16*EI_1);
        GII     = G_II_f+G_II_b;
        G       = GI+GII;
        mr      = GII./G;
    case 'Brussat1977'
        % Source:
        %   > "Stress Analysis of the Cracked Lap Shear Specimen: An ASTM
        %   round Robin" (1986) by W. Johson
        %
        % Description:
        %   > SERR for an infinitely long joint (independend of crack
        %   length
        
        G   = P.^2/(4*E*t);
        GI  = 2*M_k.^2/(7*EI_1)*(1-EI_1/EI_0);
        GII = G-GI;
    case 'Fern1und1994'
        %   > NOT FINISHED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % Source:
        %   > "Fracture load predictions for adhesive joints" (1994) by
        %   Fernlund et al.
        %
        % Description:
        %   > SERR for the J-integral and mode partitioning by Williams or
        %   Suo and Hutchinson

        G   = P.^2/(2*EA_1)+M_k.^2/(2*EI_1)-(P.^2/(2*(EA_1+EA_1))+M_0.^2/(2*EI_0));     
end

% Set output
serr.G      = G;
serr.GI     = GI;
serr.GII    = GII;
serr.MR     = serr.GII./serr.G;

end






