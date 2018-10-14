function [dbdN, dG1eq] = adhesiveDGR(GI, GII, c0, m0, c100, dbdNth)

% Mode Ratio
MR = GII./(GI+GII);

% Mode I equivalent SERR compenent
G1eq   = sqrt(GI)/2+sqrt(GI/2+GII);
dG1eq  = (G1eq(:,:,2)-G1eq(:,:,1)).^2;

% Paris Law 'C' coefficient based on the current Mode Ratio
C_mr = c100.^MR(:,:,2).*c0.^(1-MR(:,:,2));

% Get f(G) value that corresponds to the fatigue threshold
% Note: for this application, as there is rarely fatigue treshold data, the
%       user can enter a disbond growth rate that is correlated to the 
%       fatigue threshold (the vertical assymptote in the db/dN vs f(G) 
%       curve
dG1_eq_th = (dbdNth./C_mr).^(1/m0);

% Set to 0 if equal to or below the fatigue treshold
dG1eq(dG1eq<=dG1_eq_th) = 0;

% Disbond growth rate using the MR at the maximum load
dbdN = c100.^MR(:,:,2).*c0.^(1-MR(:,:,2)).*dG1eq.^m0;

% Negative dbdN is not allowed
dbdN(dbdN<0) = 0;

end