function [Minor, dMinor, dN, N, Nf] = metalFatigueInitiation(method, Sa, Sm, dbdN, dlb)

% Load-ratio
R_nom = (Sm-Sa)./(Sm+Sa);

% Pre-allocate memory
Nf = inf(size(Sa));

switch method 
    case 'MH - Sheet'
        % Source:
        %   > Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
        %   > Page 3-115; Aluminium 2024; Bare sheet, 0.090-inch (2.286 mm) thick
        S_max = Sa+Sm;
        % Equivalent maximum stress
        S_eq = S_max.*(1-R_nom).^0.56;
        % ksi to Pa
        S_eq = S_eq*1.45038e-7;
        % Cycles untill inititation (= 0 if below the threshold)
        a           = S_eq-15.8;
        Nf(a>0)    = 10.^(11.1-3.97*log10(a(a>0)));
        
    case 'MH - Rod)'
        % Source:
        %   > Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
        %   > Page 3-111; Aluminium 2024; Rolled bar
        S_max = Sa+Sm;
        % Equivalent maximum stress
        S_eq = S_max.*(1-R_nom_ad1).^0.52;
        % ksi to Pa
        S_eq = S_eq*1.45038e-7;
        % Cycles untill inititation (= 0 if below the threshold)
        Nf = 10.^(20.83-9.09*log10(S_eq));
end

%% Fatigue damage accumulation

% Cycle increment (=inf when dbdN == 0)
dNrw = dlb./dbdN;
dN = floor(dNrw);
 
% Minor Rule
dMinor = dN./Nf;

% The following scenario's and corresponding values in 'dMinor' indicate;
%   1) value = Max stress > S-N threshold   &&   db/dN > 0
%   2) 0     = Max stress < S-N threshold   &&   db/dN > 0
%   3) inf   = Max stress > S-N threshold   &&   db/dN = 0
%   4) NaN   = Max stress < S-N treshold    &&   db/dN = 0
%  Note: the fourth scenario implies infinite fatigue life...

% Accumulated total fatigue damage
Minor = nancumsum(dMinor,1,2);

% Check for scenario (3) and (4) and act accordingly
for i = 1:size(dMinor,1)
    if any(any(isinf(dMinor(i,:,:))))
        % Scenario: arrested adhesive crack
        %   3) inf   : Max stress > S-N threshold   &&   db/dN = 0

        if i == 1
            % Remaining cycles until fatigue intiation
            cycles = min(min(Nf(i,:,:)));
            dN(i)               = cycles;
            dN(i+1:end)         = 0;
            dMinor(i,:,:)       = cycles./Nf(i,:,:);
            dMinor(i+1:end,:,:) = 0;
            Minor               = nancumsum(dMinor,1,2);
        else
            if any(any(Minor(i-1,:,:)>=1))
                % Fatigue already has been initiated; no remaining cycles
                dN(i:end)           = 0;
                dMinor(i:end,:,:)   = 0;
                Minor               = nancumsum(dMinor,1,2);
            else
                % Remaining cycles until fatigue intiation
                cycles              = min(min((1-Minor(i-1,:,:)).*Nf(i,:,:)));
                dN(i)               = cycles;
                dN(i+1:end)         = 0;   
                dMinor(i,:,:)       = cycles./Nf(i,:,:);
                dMinor(i+1:end,:,:) = 0;   
                Minor               = nancumsum(dMinor,1,2);
            end
        end

        % Terminate the for loop
        break
        
    elseif any(any(isnan(Minor(i,:,:)))) && ~any(any(isinf(Minor(i,:,:))))
        % Scenario: arrested adhesive crack and no aluminum fatigue
        %           intitiation; infinite life
        %   4) NaN   : Max stress < S-N treshold    &&   db/dN = 0
        
        % Set to 2 to mark infinite fatigue life
        dN(i:end) = 2;
        Nf(i:i:end) = inf;
        dMinor(i:i:end) = 2;
        
        % Terminate the for loop
        break
    end
end

% Acumulated load cycle count
N = cumsum(dN);

end
