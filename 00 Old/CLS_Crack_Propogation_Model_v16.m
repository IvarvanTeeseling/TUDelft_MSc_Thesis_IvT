%% Matlab reset

clear;
close all;
clc;

%% User input

% Select adherent material
AdherentSelect = {'Mono_1' 'Mono_2' 'FML_1'};
AdherentSelect = AdherentSelect{3};

switch AdherentSelect
    case 'Mono_1'
        % Input from E. Verreman
        E   = 72e9;
        t   = 0.0032; % [m]
        v   = 0.33; % [-]
        L_1 = 0.0032/0.005;
        L_2 = 0.0032/0.005*0.75;
        b0  = 0;
    case 'Mono_2'
        % Monolithic adherent properties
        E   = 72e9;
        t   = 0.0035; % [m]
        v   = 0.2956; % [-]
        L_1 = 0.85;
        L_2 = 0.75;
        b0  = 0;
    case 'FML_1'
        % Dimensions from: Bachelorthesis K. Hoidus P35
        L_1 = 0.05+0.145; % Strap length (l_A + l_B)
        L_2 = 0.145; % Lap length (l_A)
        d   = 0.05; % Specimen width
        b0  = 0;
        
        % Aluminum ply properties
        Al.E1   = 72000e6;
        Al.E2   = 72000e6;
        Al.G    = 27068e6;
        Al.v12  = 0.33;
        Al.t    = 0.3e-3;
        Al.Su   = 469e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
        Al.y1   = 347e6;
        Al.y2   = 299e6;
        Al.ct1  = 2.32e-5;
        Al.ct2  = 2.32e-5;
        
        % Glass Fibre ply properties
        GF.E1   = 48900e6;
        GF.E2   = 5500e6;
        GF.G    = 5500e6;
        GF.v12  = 0.33;
        GF.t    = 0.133e-3;
        GF.ct1  = 6.1e-6;
        GF.ct2  = 26.2e-6;
        
        % Important: Must be (1) balanced, (2) symmetric and (3) quasi-isotropic
        % 1 = AL, 2 = GF
        fml.layer   = [1 2 1 2 1 2 1 2 1];
        fml.theta   = [0 0 0 90 0 90 0 0 0];
        fml.t       = (fml.layer == 1)*Al.t + (fml.layer == 2)*GF.t;
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

% Load case from: 75% from the load in Bachelorthesis K. Hoidus P49
Smax   = (12.74+10.42)*1e3/(d*t)*0.5; % [N/m^2]
Rload  = 0.1;

% Elements for crack evaluation
q       = 100;
qmax    = ceil(q*9/10);

% Support boundary conditions - SERR
BC_serr = {'RollerRoller' 'ClampedClamped'};
BC_serr = BC_serr{2};

% Discretize option
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{2};

% Select laminate stiffness option
Modulus = {'Youngs' 'Bending'};
Modulus = Modulus{1};

%% Geometric parameters

if exist('Al','var') && exist ('GF','var') && exist ('fml','var')
    % Generate the ABD matrix, Compliance matrix and laminate Young's Modulus
    % (based on (m) membrane and (b) bending) for a symmetric,
    % balanced laminate where D_12 = D_26 = 0
    [ABD, FMLm, FMLb] = ABD_Matrix_Generator(Al, GF, fml);
    
    % Laminate Modulus
    switch Modulus
        case 'Youngs'
            % Young's Modulus
            Ex  = FMLm.Ex;
            Ey  = FMLm.Ey;
            vxy = FMLm.vxy;
            vyx = FMLm.vyx;
            Gxy = FMLm.Gxy;
            E   = Ex;
            v   = vxy;
        case 'Bending'
            % Bending Modulus
            Ex  = FMLb.E1;
            Ey  = FMLb.E2;
            vxy = FMLb.v12;
            vyx = FMLb.v21;
            Gxy = FMLb.G12;
            E   = Ex;
            v   = vxy;
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
P_max   = Smax*t;               % [N/m]
P_min   = P_max*Rload;          % [N/m]
P       = [P_min P_max];        % [N/m]
P       = permute(P,[3,1,2]);   % size [1x1xn] where n = # loads

% Applied stress
S_min   = Smax*Rload;           % [N/m^2]
S       = [S_min Smax];         % [N/m^2]
S       = permute(S,[3,1,2]);   % size [1x1xn] where n = # loads

%% Discretize adherents

switch DiscretizeMethod
    case 'LeftBoundary'
        x1 = 0:dla:l_A0-dla;
        x0 = 0:dlb:l_B0-dlb;
    case 'Central'
        x1 = dla/2:dla:l_A0-dla/2;
        x0 = dlb/2:dlb:l_B0-dlb/2;
    case 'RightBoundary'
        x1 = dla:dla:l_A0;
        x0 = dlb:dlb:l_B0;
end

% Create x0 vector for each time an entire element has been cracked
x1 = x1.*ones(qmax,1);
x1 = x1 - dla*ones(qmax,q).*(0:qmax-1)';
x1 = x1.*triu(ones(size(x1)));

x0 = x0.*ones(qmax,1);
x0 = x0 - dlb*ones(qmax,q).*(0:qmax-1)';
x0 = x0.*triu(ones(size(x0)));

% Free adherent (l_A) and overlap (l_B) length
l_A = l_A0*ones(qmax,1) + dlb*(0:qmax-1)';
l_B = l_B0*ones(qmax,1) - dlb*(0:qmax-1)';

%% Overlap edge loads

% Overlap edge loads (minumum and maximum)
[M, Q] = Overlap_Edge_Loads(x1, x0, P, EIxx1, EIxx0, l_A, l_B, t, ta, BC_serr);

% Extract solutions
M_k  = M.A;
Q_k  = Q.A;
% Overlap shear force based on equilibrium
Q_kc = Q.C;

M_0 = M.B;
Q_0 = Q.B;

%% Adhesive Stresses

% x0 vector must be adjusted to the x-axis system used in the adhesive
% stress analysis: -l_B <= x0 <= 0
x00 = x0-l_B;
x00 = x00.*triu(ones(size(x00)));
%x00 = [x00 abs(fliplr(x00))];

% Adhesive stresses derrived using the horizontal force component
F = P*cos(alpha);

% Luo and Tong (2004, 2007) > adhesive thickness inlcuded
[Shear_a, Peel_a, Loads_ad1, Loads_ad2] = Calc_OverlapStressesAndLoads(x00, F, M_k, Q_kc, EIxx1, AExx1, l_B, E, t, Ea, Ga, ta);

% Goland and Reissner (1944) > adhesive thickness excluded (ta = 0)
%[Shear_a2, Peel_a2] = Adhesive_Stresses_GR(x00, P, l_B, E, t, Ea, Ga, ta, v);

%% Adherent Fatigue Accumulation

if exist('Al','var') && exist ('GF','var') && exist ('fml','var')
    % FML adherent
    
    % TO DO:
    %   1. Check error introduced by using Eb or Em (Young's or Bending
    %   Laminate Modulus)
    %   2. Check plane stress / plane strain assumptions and make sure they
    %   are included correctly in the equations
    
    % Stress and strain cycles
    [Sxx, exx] = Calc_StressStrainCycle(Loads_ad1.N, Loads_ad1.M, Loads_ad2.N, Loads_ad2.M, ABD, AExx1, EIxx1);
    
    % R-ratio (Smin/Smax)
    R_nom_ad1 = Sxx.ad1(:,:,1)./Sxx.ad1(:,:,2);
    R_nom_ad2 = Sxx.ad2(:,:,1)./Sxx.ad2(:,:,2);
    
    R_nom_ad1(isnan(R_nom_ad1)) = 0;
    R_nom_ad2(isnan(R_nom_ad2)) = 0;
    
    % Find the amplitude and mean stress of the stress cycle
    Sm_nom_ad1 = (1+R_nom_ad1)./2.*Sxx.ad1(:,:,2);
    Sa_nom_ad1 = (1-R_nom_ad1)./2.*Sxx.ad1(:,:,2);
    Sm_nom_ad2 = (1+R_nom_ad2)./2.*Sxx.ad2(:,:,2);
    Sa_nom_ad2 = (1-R_nom_ad2)./2.*Sxx.ad2(:,:,2);
    
    % Transform amplitude to a Sm = 0 cycle using the Goodman Relation
    Sa_nom_sm0_ad1 = Sa_nom_ad1./(1+Sm_nom_ad1/Al.Su);
    Sa_nom_sm0_ad2 = Sa_nom_ad2./(1+Sm_nom_ad2/Al.Su);
    
    % Temp value for now
    R_SN = 0.1;
    
    % Find equivalent S-N amplitude stress cycle
    Sa_SN_ad1 = Sa_nom_sm0_ad1./(1+Sa_nom_sm0_ad1*((2/(1-R_SN)-1)/Al.Su));
    Sa_SN_ad2 = Sa_nom_sm0_ad2./(1+Sa_nom_sm0_ad2*((2/(1-R_SN)-1)/Al.Su));
end

%% Strain Energy Release Rate

% SERR and Mode Ratio components - Specimen with finite length
[serr, mr] = Calc_StrainEnergyReleaseRate(F, M_k, M_0, E, t, AExx1, EIxx1, 'Fern1und1991');

% SERR and Mode Ratio componentes - Specimen with infinite length
[serr2, mr2] = Calc_StrainEnergyReleaseRate(F, M_k, M_0, E, t, AExx1, EIxx1, 'Brussat');

% Extract solutions
G   = serr.G;
GI  = serr.GI;
GII = serr.GII;
MR  = mr;

% Enforce DAF effect
q1      = round(1/4*qmax);
q2      = round(3/4*qmax);
daf_I   = 0.97;
daf_II  = 0.85;
x_daf   = 0:1:(q2-q1);

DAF_I   = daf_I*sin(pi/x_daf(end)*x_daf)';
DAF_II  = daf_II*sin(pi/x_daf(end)*x_daf)';

GI(q1:q2,:,:) = GI(q1:q2,:,:)-DAF_I.*GI(q1:q2,:,:);
GII(q1:q2,:,:) = GII(q1:q2,:,:)-DAF_II.*GII(q1:q2,:,:);
G = GI + GII;

MR = GII./G;

% figure(1)
% hold on
% plot(x00(1,1:qmax),GI(:,:,2));
% plot(x00(1,1:qmax),GII(:,:,2));
% plot(x00(1,1:qmax),G(:,:,2));
% hold off


%% Disbond growth rate

% Crack Growth Rate (model from D. Burger (2005), FM94 adhesive)
[dbdN] = calc_CrackGrowthRate(GI, GII, MR, 5.27*10^(-17), 3.78, 10^(-17.6));

%% Number of cycles per element

% Load cycles for the disbond to grow the element distance
n = dlb./dbdN;
% Round and store removed decimals
nrnd = floor(n);
ncut = n - nrnd;

%% Fatigue accumulation

% From 
S_eq_ad1 = Sa_nom_ad1.*(1-R_nom_ad1).^0.52;
S_eq_ad1 = S_eq_ad1*1.45038e-7;
Nf_ad1 = 10.^(20.83-9.09*log10(S_eq_ad1));

Minor = n./Nf_ad1;
MinorSum = cumsum(Minor,1);

%% Plot results

setappdata(0,'N1',Loads_ad1.N);
setappdata(0,'Q1',Loads_ad1.Q);
setappdata(0,'M1',Loads_ad1.M);
setappdata(0,'N2',Loads_ad2.N);
setappdata(0,'Q2',Loads_ad2.Q);
setappdata(0,'M2',Loads_ad2.M);
setappdata(0,'Shear_a',Shear_a);
setappdata(0,'Peel_a',Peel_a);
setappdata(0,'N',cumsum(n));
setappdata(0,'dbdN',dbdN);
setappdata(0,'MinorSum',MinorSum);
setappdata(0,'Sm_nom_ad1',Sm_nom_ad1);
setappdata(0,'Sa_nom_ad1',Sa_nom_ad1);
setappdata(0,'x',x00*1000);
setappdata(0,'Sy',Al.y1);

run('GUI_FinalPlots');

run('GUI_FinalPlots2');


% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x00(1,:),Minor(1,:))
% hold off
%
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x00(1,:),MinorSum(end,:))
% hold off

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
% plot(x00(1,:),exx.ad1(1,:,2),'r')
% plot(x00(1,:),exx.ad2(1,:,2),'b')
% hold off
% title('Strain in the xx direction')
% legend('Top adherent (bottom ply)', 'Bottom adherent (top ply)')
%
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x00(1,:),Sxx.ad1(1,:,2),'r')
% plot(x00(1,:),Sxx.ad2(1,:,2),'b')
% hold off
% title('Stress in the xx direction')
% legend('Top adherent (bottom ply)', 'Bottom adherent (top ply)')
%
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% plot(x00(end/2,:),Sxx.ad1(end/2,:,2),'r')
% plot(x00(end/2,:),Sxx.ad2(end/2,:,2),'b')
% hold off
% title('Stress in the xx direction')
% legend('Top adherent (bottom ply)', 'Bottom adherent (top ply)')


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
% title('Normal Force')
% hold on
% plot(x00(1,:),Ad_1.N(1,:,1),'r')
% plot(x00(1,:),Ad_1.N(1,:,2),'--r')
% plot(x00(1,:),Ad_2.N(1,:,1),'b')
% plot(x00(1,:),Ad_2.N(1,:,2),'--b')
% hold off
%
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% %plot(x0(1,:),Shear_a(1,:,1),'r')
% plot(x00(1,:),Shear_a(1,:,2),'b')
% title('Adhesive shear stress')
% hold off
% %
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% hold on
% %plot(x0(1,:),Peel_a(1,:,1),'r')
% plot(x00(1,:),Peel_a(1,:,2),'b')
% title('Adhesive peel stress')
% hold off
%
% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% title('Shear Force [N/m]')
% hold on
% plot(x00(1,:),Ad_1.Q(1,:,2),'r')
% plot(x00(1,:),Ad_2.Q(1,:,2),'b')
% hold off
%

% r =  findobj('type','figure');
% r = length(r);
% figure(r+1)
% title('Bending Moment [Nm^2]')
% hold on
% plot(x00(1,:),Ad_1.M(1,:,2),'r')
% plot(x00(1,:),Ad_2.M(1,:,2),'b')
% hold off
% legend('Top adherent', 'Bottom adherent')
%
%

r =  findobj('type','figure');
r = length(r);
figure(r+1)
hold on
plot(x00(1,1:qmax),GI(:,:,2),'r') % This solution
plot(x00(1,1:qmax),GII(:,:,2),'b')
plot(x00(1,1:qmax),G(:,:,2),'c')
plot(x00(1,1:qmax),serr2.G(:,:,2)*ones(1,length(x00(1,1:qmax))),'--c')
hold off

