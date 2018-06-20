function [S_xx, e_xx, R_nom, Sm_nom, Sa_nom] = Stress_Strain_Cycle(N1, M1, N2, M2, ABD, AExx, EIxx)

%% Membrane strain and curvature
%   >Note: AExx and EIxx used instead of the compliance matrix to include
%   the plane strain condition in the y-direction (if active)

% Upper adherent
e0x_1 = N1/AExx;
kx_1  = M1/EIxx;

% Lower adherent
e0x_2 = N2/AExx;
kx_2  = M2/EIxx;

% Mid-plane z-coordinate of the metal plies adjecent to the adhesive
% interface (bottom ply for top adherent, top ply for bottom adherent)
z       = ABD.zply;
z_mid   = z(1:end-1)-(z(1:end-1)-z(2:end))/2;
z1      = z_mid(1);
z2      = z_mid(end);

%% Aluminum facesheet total strain
e_xx.ad1 = e0x_1-z1*kx_1;
e_xx.ad2 = e0x_2-z2*kx_2;

%% Aluminum facesheet total stress
% TO DO: include laminate thermal stresses
% DO DO: include the notch effect
S_xx.ad1 = ABD.stiff(1,1,1)*e_xx.ad1;
S_xx.ad2 = ABD.stiff(1,1,end)*e_xx.ad2;

%% R-ratio
R_nom.ad1 = S_xx.ad1(:,:,1)./S_xx.ad1(:,:,2);
R_nom.ad2 = S_xx.ad2(:,:,1)./S_xx.ad2(:,:,2);
R_nom.ad1(isnan(R_nom.ad1)) = 0;
R_nom.ad2(isnan(R_nom.ad2)) = 0;

%% Amplitude and mean stress
Sm_nom.ad1 = (1+R_nom.ad1)./2.*S_xx.ad1(:,:,2);
Sa_nom.ad1 = (1-R_nom.ad1)./2.*S_xx.ad1(:,:,2);
Sm_nom.ad2 = (1+R_nom.ad2)./2.*S_xx.ad2(:,:,2);
Sa_nom.ad2 = (1-R_nom.ad2)./2.*S_xx.ad2(:,:,2);

end