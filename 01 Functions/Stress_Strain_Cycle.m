function [S_xx, e_xx, R_nom, Sm_nom, Sa_nom] = Stress_Strain_Cycle(N_1, M_1, N_2, M_2, ABD, AE_xx, EI_xx)
% TO DO: include laminate thermal stresses
% DO DO: include the notch effect

%% Membrane strain and curvature
%   >Note: AExx and EIxx used instead of the compliance matrix to include
%   the plane strain condition in the y-direction (if active)

% Upper adherent
e_0x_1 = N_1/AE_xx;
kx_1  = M_1/EI_xx;

% Lower adherent
e_0x_2 = N_2/AE_xx;
kx_2  = M_2/EI_xx;

% Mid-plane z-coordinate of the metal plies adjecent to the adhesive
% interface (bottom ply for top adherent, top ply for bottom adherent)
z       = ABD.zply;
z_mid   = z(1:end-1)-(z(1:end-1)-z(2:end))/2;

% Pre-allocate memory
e_xx.ad1    = zeros([size(N_1) 2]);
e_xx.ad2    = zeros([size(N_1) 2]);
S_xx.ad1    = zeros([size(N_1) 2]);
S_xx.ad2    = zeros([size(N_1) 2]);
R_nom.ad2   = zeros(size(N_1,1), size(N_1,2), 2);
R_nom.ad2   = zeros(size(N_1,1), size(N_1,2), 2);
Sm_nom.ad2  = zeros(size(N_1,1), size(N_1,2), 2);
Sm_nom.ad2  = zeros(size(N_1,1), size(N_1,2), 2);
Sa_nom.ad2  = zeros(size(N_1,1), size(N_1,2), 2);
Sa_nom.ad2  = zeros(size(N_1,1), size(N_1,2), 2);

cnt = 1;

% Stress/strain per ply
for i = [1 length(z_mid)]
    % note: e_xx.ad1 is a [NxMxPxQ] matrix where:
    %   N = the number of cracked elements
    %   M = the number of elements
    %   P = 2: the min and max load
    %   Q = 2: the two outer plies
    
    % Strain
    e_xx.ad1(:,:,:,cnt) = e_0x_1-z_mid(i)*kx_1;
    e_xx.ad2(:,:,:,cnt) = e_0x_2-z_mid(i)*kx_2;
    
    % Stress
    S_xx.ad1(:,:,:,cnt) = ABD.stiff(1,1,i)*e_xx.ad1(:,:,:,cnt);
    S_xx.ad2(:,:,:,cnt) = ABD.stiff(1,1,i)*e_xx.ad2(:,:,:,cnt);
    
    % R-ratio
    R_nom.ad1(:,:,cnt)        = S_xx.ad1(:,:,1,cnt)./S_xx.ad1(:,:,2,cnt);
    R_nom.ad2(:,:,cnt)        = S_xx.ad2(:,:,1,cnt)./S_xx.ad2(:,:,2,cnt);
    R_nom.ad1(isnan(R_nom.ad1)) = 0;
    R_nom.ad2(isnan(R_nom.ad2)) = 0;
    
    % Amplitude and mean stress
    Sm_nom.ad1(:,:,cnt) = (1+R_nom.ad1(:,:,cnt))./2.*S_xx.ad1(:,:,2,cnt);
    Sa_nom.ad1(:,:,cnt) = (1-R_nom.ad1(:,:,cnt))./2.*S_xx.ad1(:,:,2,cnt);
    Sm_nom.ad2(:,:,cnt) = (1+R_nom.ad2(:,:,cnt))./2.*S_xx.ad2(:,:,2,cnt);
    Sa_nom.ad2(:,:,cnt) = (1-R_nom.ad2(:,:,cnt))./2.*S_xx.ad2(:,:,2,cnt);
   
    cnt = cnt + 1;
    
end

end