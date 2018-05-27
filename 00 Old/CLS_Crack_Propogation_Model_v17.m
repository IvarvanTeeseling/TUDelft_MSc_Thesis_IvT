%% Matlab reset

clear;
close all;
clc;

%% Module 0: Input

% Adherent material
AdherentSelect = {'Mono_1' 'Mono_2' 'FML_1'};
AdherentSelect = AdherentSelect{3};
switch AdherentSelect
    case 'Mono_1'
        % Input from E. Verreman
        E   = 72e9;
        t   = 0.0032; % [m]
        v   = 0.33; % [-]
    case 'Mono_2'
        % Monolithic adherent properties
        E   = 72e9;
        t   = 0.0035; % [m]
        v   = 0.2956; % [-]
    case 'FML_1'
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
        
        % Lay-up (1 = AL, 2 = GF)
        fml.layer   = [1 2 1 2 1 2 1 2 1];
        fml.theta   = [0 0 0 90 0 90 0 0 0];
        fml.t       = (fml.layer == 1)*Al.t + (fml.layer == 2)*GF.t;
        t = sum(fml.t);
        
        % Laminate stiffness modulus
        Modulus = {'Youngs' 'Bending'};
        Modulus = Modulus{1};
end

% Adhesive material
AdhesiveSelect = {'FM94' 'adh_1' 'adh_2'};
AdhesiveSelect = AdhesiveSelect{3};
switch AdhesiveSelect
    case 'FM94'
        t_a      = 0.125*10^(-3);
        E_a      = 3.1e9; % (not from FM94)
        G_a      = 823e6;
        v_a      = 0.4; % (not from FM94)
        c_0      = 5.27*10^(-17);
        m_0      = 3.78;
        c_100    = 10^(-17.6);
    case 'adh_1'
        % Input from E. Verreman
        t_a  = 0.3e-3;  % [m]
        E_a  = 3.1e9; % [Pa]
        G_a  = 1.1e9; % [Pa]
        v_a  = 0.4;
        % TEMP > not real data
        c_0      = 5.27*10^(-17);
        m_0      = 3.78;
        c_100    = 10^(-17.6);
    case 'adh_2'
        t_a  = 0.2e-3;  % [m]
        E_a  = 0.7e9; % [Pa]
        G_a  = 0.7e9; % [Pa]
        v_a  = 0.4;
        % TEMP > not real data
        c_0      = 5.27*10^(-17);
        m_0      = 3.78;
        c_100    = 10^(-17.6);
end

% CLS configuration
ConfigSelect = {'Thesis' 'Config_1' 'Config_2'};
ConfigSelect = ConfigSelect{1};
switch ConfigSelect
    case 'Thesis'
        % Dimensions from Bachelorthesis K. Hoidus, page 35
        L_1 = 0.05+0.145;
        L_2 = 0.145;
        d   = 0.05;
        b_0  = 0;
    case 'Config_1'
        L_1 = 0.0032/0.005;
        L_2 = 0.0032/0.005*0.75;
        d   = 0.05;
        b_0  = 0;
    case 'Config_2'
        L_1 = 0.85;
        L_2 = 0.75;
        d   = 0.05;
        b_0  = 0;
end

% Applied load
LCSelect = {'Thesis' 'LC_1' 'LC_2'};
LCSelect = LCSelect{1};
switch LCSelect
    case 'Thesis'
        % Load case from: 75% from the load in Bachelorthesis K. Hoidus P49
        S_max   = (12.74+10.42)*1e3/(d*t)*0.5; % [N/m^2]
        R_load  = 0.1;
    case 'LC_1'
    case 'LC_2'
end

% Support boundary conditions (only for overlap loads)
BC = {'RollerRoller' 'ClampedClamped'};
BC = BC{1};

% Numerical settings - Number of elements
q       = 300;
% Numerical settings - Number of cracked elements
q_max   = ceil(q*9/10);

% Numerical settings - Discretization method
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{2};

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
AExx0 = E*t*2 + E_a*t_a;

% Bending stifness
EIxx1 = E*t^3/12;
EIxx0 = 2*E*(t^3/12+t*(t/2+t_a/2)^2)+E_a*t_a^3/12;

% Initital geometry
l_A0 = L_1-(L_2-b_0);
l_B0 = L_2-b_0;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t + t_a)./(2*(l_A0+l_B0));

% Element length
dla = l_A0/q;
dlb = l_B0/q;

%% Load parameters

% Applied force
P_max   = S_max*t;               % [N/m]
P_min   = P_max*R_load;          % [N/m]
P       = [P_min P_max];        % [N/m]
P       = permute(P,[3,1,2]);   % size [1x1xn] where n = # loads

% Applied stress
S_min   = S_max*R_load;           % [N/m^2]
S       = [S_min S_max];         % [N/m^2]
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
x1 = x1.*ones(q_max,1);
x1 = x1 - dla*ones(q_max,q).*(0:q_max-1)';
x1 = x1.*triu(ones(size(x1)));

x0 = x0.*ones(q_max,1);
x0 = x0 - dlb*ones(q_max,q).*(0:q_max-1)';
x0 = x0.*triu(ones(size(x0)));

% Free adherent (l_A) and overlap (l_B) length
l_A = l_A0*ones(q_max,1) + dlb*(0:q_max-1)';
l_B = l_B0*ones(q_max,1) - dlb*(0:q_max-1)';

%% Overlap edge loads

% Overlap edge loads (minumum and maximum)
[M, Q] = Overlap_Edge_Loads(x1, x0, P, EIxx1, EIxx0, l_A, l_B, t, t_a, BC);

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
[Shear_a, Peel_a, Loads_ad1, Loads_ad2] = Calc_OverlapStressesAndLoads(x00, F, M_k, Q_kc, EIxx1, AExx1, l_B, E, t, E_a, G_a, t_a);

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
q1      = round(1/4*q_max);
q2      = round(3/4*q_max);
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
[dbdN] = calc_CrackGrowthRate(GI, GII, MR, c_0, m_0, c_100);

%% Number of cycles per element

% Load cycles for the disbond to grow the element distance
n = dlb./dbdN;
% Round and store removed decimals
nrnd = floor(n);
ncut = n - nrnd;

%% Fatigue accumulation

% Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
% Page 115; Aluminium 2024; Bare sheet, 0.090-inch (2.286 mm) thick
% S_eq_ad1 = Sa_nom_ad1.*(1-R_nom_ad1).^0.56;     % Equivalent amplitude S
% S_eq_ad1 = S_eq_ad1*1.45038e-7;                 % Pa to ksi
% Nf_ad1 = 10.^(11.1-3.97*log10(S_eq_ad1-15.8));  % Nr. cycles untill fatigue initiation

% From 
S_eq_ad1 = Sa_nom_ad1.*(1-R_nom_ad1).^0.52;
S_eq_ad1 = S_eq_ad1*1.45038e-7;
Nf_ad1 = 10.^(20.83-9.09*log10(S_eq_ad1));

% Fatigue damage accumulation
Minor = n./Nf_ad1;                              % Damage per cracked element
MinorSum = cumsum(Minor,1);                     % Accumulated total damage

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
plot(x00(1,1:q_max),GI(:,:,2),'r') % This solution
plot(x00(1,1:q_max),GII(:,:,2),'b')
plot(x00(1,1:q_max),G(:,:,2),'c')
plot(x00(1,1:q_max),serr2.G(:,:,2)*ones(1,length(x00(1,1:q_max))),'--c')
hold off

