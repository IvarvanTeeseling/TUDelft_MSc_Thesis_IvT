%% Matlab reset

clear;
close all;
clc;

%% User input

% Select adherent material
AdherentSelect = {'Mono_1' 'Mono_2' 'FML_1' 'FML_2'};
AdherentSelect = AdherentSelect{2};

switch AdherentSelect
    
    case 'Mono_1'
        % Input from E. Verreman
        E   = 72e9;
        t   = 0.0032; % [m]
        v   = 0.33; % [-]
        L_1 = 0.25;
        L_2 = 0.2;
        b0  = 0;
    case 'Mono_2'
        % Monolithic adherent properties
        E   = 70e9;
        t   = 0.004; % [m]
        v   = 0.34; % [-]
        L_1 = 0.05;
        L_2 = 0.01;
        b0  = 0;
    case 'FML_1'
        L_1 = 0.15;
        L_2 = 0.1;
        b0  = 0;
        
        % Aluminum ply properties
        Al.Ex   = 72400e6;
        Al.Ey   = 72400e6;
        Al.G    = 26700e6;
        Al.vxy  = 0.33;
        Al.t    = 0.3e-3;
        
        % Glass Fibre ply properties
        GF.Ex   = 48900e6;
        GF.Ey   = 5500e6;
        GF.G    = 5550;
        GF.vxy  = 0.33;
        GF.t    = 0.133*10^(-3);
        
        % Important: Must be (1) balanced, (2) symmetric and (3) quasi-isotropic
        fml.layer   = [1 2 1 2 2 1 2 1 1 2 1 2 2 1 2 1];
        fml.t       = [Al.t GF.t Al.t GF.t GF.t Al.t GF.t Al.t Al.t GF.t Al.t GF.t GF.t Al.t GF.t Al.t];
        fml.theta   = [0 90 0 0 0 0 90 0 0 90 0 0 0 0 90 0];
        t = sum(fml.t); 
            case 'FML_2'
        % Aluminum ply properties
        Al.Ex   = 72400e6;
        Al.Ey   = 72400e6;
        Al.G    = 26700e6;
        Al.vxy  = 0.33;
        Al.t    = 0.3e-3;
        
        % Glass Fibre ply properties
        GF.Ex   = 48900e6;
        GF.Ey   = 5500e6;
        GF.G    = 5550;
        GF.vxy  = 0.33;
        GF.t    = 0.133*10^(-3);
        
        % Important: Must be (1) balanced, (2) symmetric and (3) quasi-isotropic
        fml.layer   = [1 2 1 2 2 1 2 1];
        fml.t       = [Al.t GF.t Al.t GF.t GF.t Al.t GF.t Al.t];
        fml.theta   = [0 90 0 0 0 0 90 0];
        t = sum(fml.t);
end

% Select adhesive material
AdhesiveSelect = {'FM94' 'adh_1' 'adh_2'};
AdhesiveSelect = AdhesiveSelect{3};

switch AdhesiveSelect
    case 'FM94'
        ta = 0.125*10^(-3);
        Ea = 3.1e9; % (not from FM94)
        Ga = 823e6;
        va = 0.4; % (not from FM94)
    case 'adh_1'
        % Input from E. Verreman
        ta  = 0.3e-3;  % [m]
        Ea  = 3.1e9; % [Pa]
        Ga  = 1.1e9; % [Pa]
        va  = 0.4;
    case 'adh_2'
        ta  = 0.2e-3;  % [m]
        Ea  = 0.7e9; % [Pa]
        Ga  = 0.7e9; % [Pa]
        va  = 0.4;
end

% Load case
Smax   = 60e3/t; % [N/m^2]
Rload  = 0.1;

% Support boundary conditions - SERR
BC_serr = {'RollerRoller' 'ClampedClamped'};
BC_serr = BC_serr{1};

% Discretize option
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{1};

% Select laminate stiffness option
Modulus = {'Youngs' 'Bending'};
Modulus = Modulus{1};

% Elements for crack evaluation
q = 50;
qmax = ceil(q*10/10);

%% Geometric parameters

if exist('Al','var') && exist ('GF','var') && exist ('fml','var')
    % Generate the ABD matrix, Compliance matrix and laminate Young's Modulus
    % (based on (m) membrane and (b) bending) for a symmetric,
    % balanced laminate where D_12 = D_26 = 0
    [ABD, FMLm, FMLb] = ABD_Matrix_Generator(Al, GF, fml);
    
    % FML thickness
    t = sum(fml.t);
    
    % Laminate Modulus
    switch Modulus
        case 'Youngs'
            % Young's Modulus
            E1  = FMLm.E1;
            E2  = FMLm.E2;
            v12 = FMLm.v12;
            v21 = FMLm.v21;
            G12 = FMLm.G12;
            % TO DO: check if E1 ~= E2 is a problem w.r.t. the solution
            E = E1;
            v = v12;
        case 'Bending'
            % Bending Modulus
            E1  = FMLb.E1;
            E2  = FMLb.E2;
            v12 = FMLb.v12;
            v21 = FMLb.v21;
            G12 = FMLb.G12;
            
            E = E1;
            v = v12;
    end
end

% Membrane stiffness
AExx1 = E*t;
AExx0 = E*t*2 + Ea*ta;

% Bending stifness
EIxx1 = E*t^3/12;
EIxx0 = 2*E*(t^3/12+t*(t/2+ta/2)^2)+Ea*ta^3/12;

% Initital geometry
l_A0 = L_1-(L_2-b0);
l_B0 = L_2-b0;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t + ta)./(2*(l_A0+l_B0));

% Element length
dla = l_A0/q;
dlb = l_B0/q;

%% Load parameters

% Applied force
P_max   = Smax*t;          % [N/m]
P_min   = P_max*Rload;     % [N/m]
P       = [P_min P_max];   % [N/m]
P       = permute(P,[3,1,2]);

% Applied stress
S_min   = Smax*Rload;     % [N/m^2]
S       = [S_min Smax];    % [N/m^2]
S       = permute(S,[3,1,2]);

%% Discretize adherents

switch DiscretizeMethod
    case 'LeftBoundary'
        x1 = 0:dla:l_A0-dla;
        x0 = 0:dlb:l_B0-dlb;
        %x0 = [x0 x0(:,end)+dlb];
    case 'Central'
        x1 = dla/2:dla:l_A0-dla/2;
        x0 = dlb/2:dlb:l_B0-dlb/2;
        %x0 = [x0 x0(:,end)+dlb/2];
    case 'RightBoundary'
        x1 = dla:dla:l_A0;
        x0 = dlb:dlb:l_B0;
end

%% Overlap edge loads

% Create x0 vector for each time an entire element has been cracked
x1 = x1.*ones(qmax,1);
x1 = x1 - dla*ones(qmax,q).*(0:qmax-1)';
x0 = x0.*ones(qmax,1);
x0 = x0 - dlb*ones(qmax,q).*(0:qmax-1)';
% Remove floating error (don't know the source)
x1(x1<1e-15) = 0;
x0(x0<1e-15) = 0;
for j = 1:size(x0,1)-1
    % All entries up to the first element must equal the first element
    % Ugly, but functional code..
    x1(j+1,1:j) = x1(j+1,j+1);
    x0(j+1,1:j) = x0(j+1,j+1);
end

% Free adherent (l_A) and overlap (l_B) length
l_A = l_A0*ones(qmax,1) + dlb*(0:qmax-1)';
l_B = l_B0*ones(qmax,1) - dlb*(0:qmax-1)';

% Overlap edge loads (minumum and maximum)
[M, Q] = Overlap_Edge_Loads(x1, x0, P, EIxx1, EIxx0, l_A, l_B, t, ta, BC_serr);

% Extract solutions
M_k = M.A;
Q_k = Q.A;

M_0 = M.B;
Q_0 = Q.B;

%% Adhesive Stresses

% x0 vector must be adjusted to the x-axis system used in the adhesive
% stress analysis: -l_B <= x0 <= 0
x00 = x0-l_B;
x00 = [x00 abs(fliplr(x00))];

% Adhesive stresses derrived using the horizontal force component
F = P*cos(alpha);

% Luo and Tong (2004, 2007) > adhesive thickness inlcuded
[Shear_a, Peel_a, Ad_1, Ad_2] = Adhesive_Stresses_LuoTong(x00, F, M_k, Q_k, EIxx1, AExx1, l_A, l_B, E, t, Ea, Ga, ta);

% Goland and Reissner (1944) > adhesive thickness excluded (ta = 0)
[Shear_a2, Peel_a2] = Adhesive_Stresses_GR(x00, P, l_B, E, t, Ea, Ga, ta, v);

%% Adherent Fatigue Accumulation

if exist('Al','var') && exist ('GF','var') && exist ('fml','var')
    % FML adherent
    
    % TO DO:
    %   1. Check error introduced by using Eb or Em (Young's or Bending
    %   Laminate Modulus)
    %   2. Check plane stress / plane strain assumptions and make sure they
    %   are included correctly in the equations
    
    % Membrane strain and curvature - Upper adherent
    e0x_1 = Ad_1.N/AExx1;
    kx_1  = Ad_1.M/EIxx1;
    % Membrane strain and curvature - Lower adherent
    e0x_2 = Ad_2.N/AExx1;
    kx_2  = Ad_2.M/EIxx1;
    
    % Membrane strain and curvature
    %e   = ABD.ABD(1:3,1:3)^(-1)*[Ad_1.N(1,1,1) 0 0]';
    %k   = ABD.ABD(4:6,4:6)^(-1)*[Ad_1.M(1,1,1) 0 0]';
    %e0x = e(1);
    %kx  = k(1);
    
    % Mid-plane of the aluminum ply adjecent to bonded line interface
    z1 = ABD.zply(1:end-1)-(ABD.zply(1:end-1)-ABD.zply(2:end))/2;
    z2 = ABD.zply(1:end-1)-(ABD.zply(1:end-1)-ABD.zply(2:end))/2;
    
    % Bottom ply for adherent 1 and top ply for adherent 2
    z1 = z1(1);
    z2 = z2(end);
    
    % Laminate strain at the ply mid plane coordinate
    exx1 = e0x_1 - z1*kx_1;
    exx2 = e0x_2 - z2*kx_2;
    
    exx1m = e0x_1;
    exx2m = e0x_2;
    exx1b = -z1*kx_1;
    exx2b = -z1*kx_2;
    
    % Total laminate stress (plane strain)
    s1 = ABD.stiff(:,:,1)*[exx1(1) 0 0]';
    s2 = ABD.stiff(:,:,1)*[exx1(2) 0 0]';
end

%% Strain Energy Release Rate

% Lai et al. (1995); using the Suo and Hutchinson (1990) approach
mu = 1;
sm = 1;
delta = 1;
I = sm*((delta-1/mu)^2-(delta-1/mu)+1/3)+delta/mu*(delta-1/mu)+1/(3*mu^3);
A = 1/mu + sm;
G_T1 = 1/(2*E)*((P.^2/t+12*M_k.^2./t^3)+(-P.^2/(A*t)-M_0.^2/(I*t^3)));

% MR of 5 assumed from FEM data
G_I1 = G_T1/6;
G_II1 = G_I1*5;
MR = G_II1./G_I1;

G_II = G_II1;
G_I = G_I1;
G_T = G_T1;

% Brussat et al (1977); infinitely long specimen 
G_T2 = P.^2/(4*E*t);

% SERR
%G_I  = ta/(2*Ea)*Peel_a(:,1,:).^2;
%G_II = ta/(2*Ga)*Shear_a(:,1,:).^2;
%G_T  = G_I+G_II;

% Mode ratio
% Note: MR different at P_min and _P_max due to geometric non linearity
MR  = G_II./G_T;

%% Disbond growth rate 
% Applying D. Burger (2005) PhD data and model for FM94 adhesive

% Equivalent mode I cycle range
G1_eq = sqrt(G_I)/2+sqrt(G_I/2+G_II);
dG1_eq = (G1_eq(:,:,2)-G1_eq(:,:,1)).^2;

% Paris Law fit coefficients (D. Burger data on FM94)
% TO DO: include threshold
% TO DO: G_c
c0 = 5.27*10^(-17);
m0 = 3.78;
c100 = 10^(-17.6);

% Disbond growth rate using the MR at the maximum load
dbdN = c100.^MR(:,:,2).*c0.^(1-MR(:,:,2)).*dG1_eq.^m0;

% Load cycles for the disbond to grow the element distance
n = dlb./dbdN;
nrnd = floor(n);
ncut = n - nrnd;

%% Fatigue accumulation



%% Plot results 


% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x00(1,:),exx1m(1,:,2),'r')
% plot(x00(1,:),exx2m(1,:,2),'b')
% plot(x00(1,:),exx1b(1,:,2),'--r')
% plot(x00(1,:),exx2b(1,:,2),'--b')
% hold off
% title('Axial strain')
% legend('Top adherent axial', 'Bottom adherent axial', 'Top adherent bending', 'Bottom adherent bending')
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x00(1,:),exx1(1,:,2),'r')
% plot(x00(1,:),exx2(1,:,2),'b')
% hold off
% title('Axial strain')
% legend('Top adherent', 'Bottom adherent')

% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% plot(cumsum(n),x0(1,1:qmax)+dlb)
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% plot(cumsum(n),dbdN)
% 
r =  findobj('type','figure');
r = length(r);
figure(r+1)
title('Normal Force')
hold on
plot(x00(1,:),Ad_1.N(1,:,1),'r')
plot(x00(1,:),Ad_1.N(1,:,2),'--r')
plot(x00(1,:),Ad_2.N(1,:,1),'b')
plot(x00(1,:),Ad_2.N(1,:,2),'--b')
hold off
% 
r =  findobj('type','figure');
r = length(r);
figure(r+1)
hold on
%plot(x0(1,:),Shear_a(1,:,1),'r')
plot(x00(1,:),Shear_a(1,:,2)/(P(1,1,2)/t),'b')
title('Adhesive shear stress')
hold off
%
r =  findobj('type','figure');
r = length(r);
figure(r+1)
hold on
%plot(x0(1,:),Peel_a(1,:,1),'r')
plot(x00(1,:),Peel_a(1,:,2)/(P(1,1,2)/t),'b')
title('Adhesive peel stress')
hold off
% 
r =  findobj('type','figure');
r = length(r);
figure(r+1)
title('Shear Force [N/m]')
hold on
plot(x00(1,:),Ad_1.Q(1,:,2),'r')
plot(x00(1,:),Ad_2.Q(1,:,2),'b')
hold off
% 
r =  findobj('type','figure');
r = length(r);
figure(r+1)
title('Bending Moment [Nm^2]')
hold on
plot(x00(1,:),Ad_1.M(1,:,1),'r')
plot(x00(1,:),Ad_2.M(1,:,1),'b')
plot(x00(1,:),Ad_1.M(1,:,2),'--r')
plot(x00(1,:),Ad_2.M(1,:,2),'--b')
hold off
legend('Top adherent min','Bottom adherent min','Top adherent max','Bottom adherent max')
% 
% 

% Data_1
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(-1*(x0(1,1:qmax)*1000-200),G_I(:,:,2),'r') % This solution
% plot(-1*(x0(1,1:qmax)*1000-200),G_II(:,:,2),'b')
% plot(-1*(x0(1,1:qmax)*1000-200),G_T(:,:,2),'c')
% plot(data1(1,:), data1(2,:), '--r') % FEM solution
% plot(data2(1,:), data2(2,:), '--b')
% plot(data3(1,:), data3(2,:), '--c')
% plot(data4(1,:), data4(2,:), ':r') % Verreman solution
% plot(data5(1,:), data5(2,:), ':b')
% plot(data6(1,:), data6(2,:), ':c')
% hold off

