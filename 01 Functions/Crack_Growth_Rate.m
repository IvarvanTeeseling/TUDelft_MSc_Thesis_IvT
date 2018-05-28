function [dbdN, dG1_eq] = Crack_Growth_Rate(GI, GII, MR, c0, m0, c100)
% TO DO: include threshold
% TO DO: G_c
% TO DO: check if the correct values

% Equivalent Mode I cycle range
G1_eq   = sqrt(GI)/2+sqrt(GI/2+GII);
dG1_eq  = (G1_eq(:,:,2)-G1_eq(:,:,1)).^2;

% Disbond growth rate using the MR at the maximum load
dbdN = c100.^MR(:,:,2).*c0.^(1-MR(:,:,2)).*dG1_eq.^m0;

end