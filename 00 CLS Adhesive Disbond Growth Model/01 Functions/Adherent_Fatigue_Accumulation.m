function [Minor_csm, Minor, dN, N_f] = Adherent_Fatigue_Accumulation(method, Sa_nom, Sm_nom, Su, dbdN, dlb)

% Load-ratio
R_nom = (Sm_nom+Sa_nom)./(Sm_nom-Sa_nom);

% Pre-allocate memory
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
        % ksi to Pa
        S_eq = S_eq*1.45038e-7;
        % Cycles untill inititation (= 0 if below the threshold)
        a           = S_eq-15.8;
        N_f(a>0)    = 10.^(11.1-3.97*log10(a(a>0)));
        
    case 'Military Handbook - Rod)'
        % Source:
        %   > Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
        %   > Page 3-111; Aluminium 2024; Rolled bar
        S_max = Sa_nom+Sm_nom;
        % Equivalent maximum stress
        S_eq = S_max.*(1-R_nom_ad1).^0.52;
        % ksi to Pa
        S_eq = S_eq*1.45038e-7;
        % Cycles untill inititation (= 0 if below the threshold)
        N_f = 10.^(20.83-9.09*log10(S_eq));
end

%% Fatigue damage accumulation

% Cycle increment (=inf when dbdN <= 0)
dN_rw           = inf(size(dbdN));
dN_rw(dbdN>0)   = dlb./dbdN(dbdN>0);
dN              = floor(dN_rw);

% Minor Rule
Minor = dN./N_f;

% The following scenario's and corresponding values in 'Minor' indicate;
%   1) value = Max stress > S-N threshold   &&   db/dN > 0
%   2) 0     = Max stress < S-N threshold   &&   db/dN > 0
%   3) inf   = Max stress > S-N threshold   &&   db/dN = 0
%   4) NaN   = Max stress < S-N treshold    &&   db/dN = 0
%  Note: the fourth implies infinite fatigue life...

% Restore 0 values for the cracked elements
N_f(tril(ones(size(N_f)),-1)==1) = 0;
Minor(tril(ones(size(Minor)),-1)==1) = 0;

% Accumulated total fatigue damage
Minor_csm = cumsum(Minor,1);

% Check for scenario (3) and (4) and act accordingly
for i = 1:size(Minor_csm,1)
    if isinf(Minor_csm(i,i))
        % Scenario: arrested adhesive crack
        %   3) inf   : Max stress > S-N threshold   &&   db/dN = 0
        
        % Remove the subsequent rows as the crack will not grow any further
        dN      = dN(1:i);
        N_f     = N_f(1:i,:);
        Minor   = Minor(1:i,:);

        % Remaining cycles until fatigue intiation
        if i == 1
            cycles = N_f(i,i);
        else
            cycles = (1-Minor_csm(i-1,i))*N_f(i,i);
        end
        
        if cycles > 0
            % Fatigue has not yet initiated
            dN(i)           = cycles;
            Minor(i,i:end)  = cycles./N_f(i,i:end);
            Minor_csm       = cumsum(Minor,1);
        else
            % Fatigue has already initiated (no remaining cycles left)
            dN(i)           = 0;
            Minor(i,i:end)  = 0;
            Minor_csm       = cumsum(Minor,1);
        end
        
        % Terminate the for loop
        break
        
    elseif isnan(Minor_csm(i,i))
        % Scenario: arrested adhesive crack and no aluminum fatigue
        %           intitiation; infinite life
        %   4) NaN   : Max stress < S-N treshold    &&   db/dN = 0
        
        % Set to 2 to mark infinite fatigue life
        dN(i:end)       = 2;
        N_f(i:i:end)    = inf;
        Minor(i:i:end)  = 2;
        
        % Terminate the for loop
        break
    end
    
end

end
