function [S_xx, e_xx, R_nom, Sm_nom, Sa_nom] = Stress_Strain_Cycle(N, M, ABD, AE_xx, EI_xx)
% TO DO: include laminate thermal stresses
% DO DO: include the notch effect

%% Membrane strain and curvature
%   >Note: AExx and EIxx used instead of the compliance matrix to include
%   the plane strain condition in the y-direction (if active)

% Membrane strain and curvature
e_x0    = N/AE_xx;
k_x     = M/EI_xx;

% Ply mid-plane coordinates
z       = ABD.zply;
z_mid   = z(1:end-1)-(z(1:end-1)-z(2:end))/2;

% Pre-allocate memory
e_xx    = zeros([size(M) 2]);
S_xx    = zeros([size(M) 2]);
R_nom   = zeros(size(M,1), size(M,2), 2);
Sm_nom  = zeros(size(M,1), size(M,2), 2);
Sa_nom  = zeros(size(M,1), size(M,2), 2);

% Initialize counter
cnt = 1;

for i = [1 length(z_mid)]
    % Note: e_xx and s_xx are a [NxMxPxQ] matrix where:
    %   N = the number of cracked elements
    %   M = the number of elements
    %   P = 2: the min (=1) and max (=2) (load
    %   Q = 2: the bottom (=1) and top (=2) ply
    
    % Total strain
    e_xx(:,:,:,cnt) = e_x0-z_mid(i)*k_x;
    
    % Total stress
    S_xx(:,:,:,cnt) = ABD.stiff(1,1,i)*e_xx(:,:,:,cnt);
    
    % R-ratio
    R_nom(:,:,cnt)      = S_xx(:,:,1,cnt)./S_xx(:,:,2,cnt);
    R_nom(isnan(R_nom)) = 0;
    
    % Amplitude and mean stress
    Sm_nom(:,:,cnt) = (1+R_nom(:,:,cnt))./2.*S_xx(:,:,2,cnt);
    Sa_nom(:,:,cnt) = (1-R_nom(:,:,cnt))./2.*S_xx(:,:,2,cnt);
   
    % Update counter
    cnt = cnt + 1;   
end

end