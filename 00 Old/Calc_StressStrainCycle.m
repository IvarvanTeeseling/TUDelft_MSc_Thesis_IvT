function [sxx, exx] = Calc_StressStrainCycle(N1, M1, N2, M2, ABD, AExx, EIxx)

% Membrane strain and curvature - Upper adherent
e0x_1 = N1/AExx;
kx_1  = M1/EIxx;
% Membrane strain and curvature - Lower adherent
e0x_2 = N2/AExx;
kx_2  = M2/EIxx;

% Mid-plane z-coordinate of the metal plies adjecent to the adhesive
% interface (bottom ply for top adherent, top ply for bottom adherent)
z1 = ABD.zply(1:end-1)-(ABD.zply(1:end-1)-ABD.zply(2:end))/2;
z2 = ABD.zply(1:end-1)-(ABD.zply(1:end-1)-ABD.zply(2:end))/2;
z1 = z1(1);
z2 = z2(end);

% Laminate strain total strain
exx.ad1 = e0x_1-z1*kx_1;
exx.ad2 = e0x_2-z2*kx_2;

% Laminate total stress
% TO DO: include laminate thermal stresses
% DO DO: include the notch effect
sxx.ad1 = ABD.stiff(1,1,1)*exx.ad1;
sxx.ad2 = ABD.stiff(1,1,end)*exx.ad2;

end