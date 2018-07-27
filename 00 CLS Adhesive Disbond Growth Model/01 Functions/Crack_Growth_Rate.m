function [dbdN, dG1_eq] = Crack_Growth_Rate(GI, GII, MR, c0, m0, c100, dadN_th)

% Equivalent Mode I cycle range
G1_eq   = sqrt(GI)/2+sqrt(GI/2+GII);
dG1_eq  = (G1_eq(:,:,2)-G1_eq(:,:,1)).^2;

% Paris Law 'C' coefficient based on the Mode Ratio
C_mr = c100.^MR(:,:,2).*c0.^(1-MR(:,:,2));

% Fatigue growth threshold
dG1_eq_th = (dadN_th./C_mr).^(1/m0);

% Set to 0 if equal to or below the fatigue treshold
dG1_eq(dG1_eq<=dG1_eq_th) = 0;

% Disbond growth rate using the MR at the maximum load
dbdN = c100.^MR(:,:,2).*c0.^(1-MR(:,:,2)).*dG1_eq.^m0;

end