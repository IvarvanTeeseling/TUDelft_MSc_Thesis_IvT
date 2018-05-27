%% Matlab reset

clear all;
close all;
clc;

%% User input

% Adherent length
%L_1     = 0.65;     % [m] % Makes an error, but do not know why
%L_2     = 0.5;      % [m]
L_1 = 0.25;         % [m]
L_2 = 0.2;          % [m]
b0  = 0;            % [m]

% Select adherent material
AdherentSelect = {'Mono_1' 'FML_1'};
AdherentSelect = AdherentSelect{1};

switch AdherentSelect
    case 'Mono_1'
        % Monolithic adherent properties
        E = 72e9;
        t = 3.2e-3;
        v = 0.33;
    case 'FML_1'
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
AdhesiveSelect = {'FM94' 'adh_1'};
AdhesiveSelect = AdhesiveSelect{1};

switch AdhesiveSelect
    case 'FM94'
        ta = 0.125*10^(-3);
        Ea = 3.1e9; % (not from FM94)
        Ga = 823e6;
        va = 0.4; % (not from FM94)
    case 'adh_1'
        ta  = 0.3e-3;
        Ea  = 3.1e9;
        Ga  = 1.1e9;
        va  = 0.4;      
end

% Load case
Smax   = 250e3/t; % [N/m^2]
Rload  = 0.1;

% Support boundary conditions
BC = {'RollerRoller' 'ClampedClamped'};
BC = BC{1};

% Discretize option
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{1};

% Select laminate stiffness option
Modulus = {'Youngs' 'Bending'};
Modulus = Modulus{2};

% Elements for crack evaluation
q = 10;
qmax = ceil(q*6/8);

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
P       = [P_min P_max];    % [N/m]

% Applied stress
S_min   = Smax*Rload;     % [N/m^2]
S       = [S_min Smax];    % [N/m^2]

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
x0 = x0.*ones(qmax,1)
x0 = x0 - dlb*ones(qmax,q).*(0:qmax-1)'
% x0(x0<0) = 0
for j = 1:size(x0,1)-1
    % All entries up to the first element must equal the first element
    % Ugly, but functional code..
    x0(j+1,1:j) = x0(j+1,j+1);
end
x0

% Free adherent (l_A) and overlap (l_B) length
l_A = l_A0*ones(qmax,1) + dlb*(0:qmax-1)';
l_B = l_B0*ones(qmax,1) - dlb*(0:qmax-1)';

% Overlap edge loads (minumum and maximum)
[M, Q] = Overlap_Edge_Loads(x0, P, EIxx1, EIxx0, l_A, l_B, t, ta, BC, 3);

% Extract solutions
M_k = M.A;
Q_k = Q.A;

M_0 = M.B;
Q_0 = Q.B;

%% Adhesive Stresses

% x0 vector must be adjusted to the x-axis system used in the adhesive
% stress analysis: -l_B <= x0 <= 0
x00 = x0-l_B;

% Adhesive stresses derrived using the horizontal force component
F = P*cos(alpha);

[Shear_a, Peel_a, Ad_1, Ad_2] = Adhesive_Stresses_LuoTong(x00, F, M_k, Q_k, EIxx1, AExx1, l_B, E, t, Ea, Ga, ta, BC);

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
    
    
    e0x_2
    z2*kx_2
    
    exx2 = e0x_2 - z2*kx_2;
    
    % Total laminate stress (plane strain)
    s1 = ABD.stiff(:,:,1)*[exx1(1) 0 0]';
    s2 = ABD.stiff(:,:,1)*[exx1(2) 0 0]';
else
    % Monolithic adherent
end

%% Strain Energy Release Rate

% Lai et al. (1995); using the Suo and Hutchinson (1990) approach
mu = 1;
sm = 1;
delta = 1;
I = sm*((delta-1/mu)^2-(delta-1/mu)+1/3)+delta/mu*(delta-1/mu)+1/(3*mu^3);
A = 1/mu + sm;
G_T1 = 1/(2*E)*((P.^2/t+12*M_k.^2./t^3)+(-P.^2/(A*t)-M_0.^2/(I*t^3)));

% Brussat et al (1977); infinitely long specimen 
G_T2 = P.^2/(4*E*t);

% SERR
G_I  = ta/(2*Ea)*Peel_a(:,1,:).^2;
G_II = ta/(2*Ga)*Shear_a(:,1,:).^2;
G_T  = G_I+G_II;

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
% plot(cumsum(n),x0(1,1:qmax)+dlb)
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% plot(cumsum(n),dbdN)
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x0(1,:),Ad_1.N(1,:,1),'r')
% plot(x0(1,:),Ad_1.N(1,:,2),'b')
% hold off
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1))
% hold on
% plot(x0(1,:),Shear_a(1,:,1),'r')
% plot(x0(1,:),Shear_a(1,:,2),'b')
% hold off
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x0(1,:),Ad_1.Q(1,:,1),'r')
% plot(x0(1,:),Ad_2.Q(1,:,1),'b')
% hold off
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x0(1,:),Ad_1.M(1,:,1),'r')
% plot(x0(1,:),Ad_2.M(1,:,1),'b')
% hold off
% 
% 
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(-1*(x0(1,1:qmax)*1000-200),G_I(:,:,2),'r') % This solution
% plot(-1*(x0(1,1:qmax)*1000-200),G_II(:,:,2),'b')
% plot(-1*(x0(1,1:qmax)*1000-200),G_T(:,:,2),'c')
% plot(-1*(x0(1,1:qmax)*1000-200),G_T1(:,2),'g')
% hold off

