function [Minor_csm, Minor, N_f] = Adherent_Fatigue_Accumulation(method, Sa_nom, Sm_nom, R_nom, Su, dN)

N_f = inf(size(Sa_nom));

switch method
    case 'S. Spronk'
        % Transform amplitude to a Sm = 0 cycle using the Goodman Relation
        Sa_nom_sm0 = Sa_nom_ad1./(1+Sm_nom_ad1/Su);
        % Temp value for now
        R_SN = 0.1;
        % Find equivalent S-N amplitude stress cycle
        Sa_SN_ad1 = Sa_nom_sm0./(1+Sa_nom_sm0*((2/(1-R_SN)-1)/Su));
        
    case 'Military Handbook - Sheet'
        % Source:
        %   > Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
        %   > Page 3-115; Aluminium 2024; Bare sheet, 0.090-inch (2.286 mm) thick
        S_max = Sa_nom+Sm_nom;
        % Equivalent maximum stress
        S_eq = S_max.*(1-R_nom).^0.56;
        % Pa to ksi
        S_eq = S_eq*1.45038e-7;
        % N_f = 0 if the maximum stress is below the threshold
        a               = S_eq-15.8;
        N_f(a>0)        = 10.^(11.1-3.97*log10(a(a>0)));
        
    case 'Military Handbook - Rod)'
        % Source:
        %   > Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
        %   > Page 3-111; Aluminium 2024; Rolled bar
        S_max = Sa_nom+Sm_nom;
        % Equivalent maximum stress 
        S_eq = S_max.*(1-R_nom_ad1).^0.52;  
        % Pa to ksi
        S_eq = S_eq*1.45038e-7;                          
        % Nr. cycles untill fatigue initiation
        N_f = 10.^(20.83-9.09*log10(S_eq));             
end

%% Fatigue damage accumulation

% The following values indicate;
%   > 0     = zero Al fatigue crack growth     & non-zero disbond growth
%   > inf   = non-zero AL fatigue crack growth & zero disbond growth
%   > NaN   = zero Al fatigue crack growth     & zero disbond growth
%  Note: the last implies infinite fatigue life...
Minor = dN./N_f;

% Restore 0 values for the cracked elements
for i = 2:size(N_f,1)
   N_f(i,1:i-1)     = 0;
   Minor(i,1:i-1)   = 0;
end

% Accumulated total damage
Minor_csm = cumsum(Minor,1);

% Check for zero Al fatigue CGR and/or zero disbond growth


end
