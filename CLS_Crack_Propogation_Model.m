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
        fml.theta   = [0 0 0 0 0 0 0 0 0];
        fml.t       = (fml.layer == 1)*Al.t + (fml.layer == 2)*GF.t;
        t = sum(fml.t);
        
        % Laminate stiffness modulus
        Modulus = {'Youngs' 'Flexural'};
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
        L_1 = 0.25+0.145;
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
LCSelect = LCSelect{2};
switch LCSelect
    case 'Thesis'
        % Load case from: 75% from the load in Bachelorthesis K. Hoidus P49
        S_max   = (12.74+10.42)*1e3/(d*t)*0.70; % [N/m^2]
        R_load  = 0.1;
    case 'LC_1'
        % Load case from R. Hanx
        S_max   = 180*1e3/t*1.5;
        R_load  = 0.1;
    case 'LC_2'
end

% Support boundary conditions (only for overlap loads)
BC = {'RollerRoller' 'ClampedClamped'};
BC = BC{2};

% Numerical settings - Number of elements
q = 1000;
% Numerical settings - Number of cracked elements
q_max = ceil(q*8/10);
% Numerical settings - Discretization method
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{2};

%% Module 1: Geometry and Load Paramaters

% Initital geometry
l_A0 = L_1-(L_2-b_0);
l_B0 = L_2-b_0;

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t+t_a)./(2*(l_A0+l_B0));

% Applied force
P_max   = S_max*t;              % [N/m]
P_min   = P_max*R_load;         % [N/m]
P       = [P_min P_max];        % [N/m]
P       = permute(P,[3,1,2]);   % size [1x1xn] where n = # loads

% Applied stress
S_min   = S_max*R_load;         % [N/m^2]
S       = [S_min S_max];        % [N/m^2]
S       = permute(S,[3,1,2]);   % size [1x1xn] where n = # loads

%% Module 2: FML Laminate Properties

if exist('Al','var') && exist ('GF','var') && exist ('fml','var')
    % Generate (1) the ABD matrix, (2) Compliance matrix and (3) laminate equivalent Modulus 
    % The laminate equivalent modulus is given for a symmetric,
    % balanced laminate where A_12 = 0, A_26 = 0, D_12 ~ 0, D_26 ~ 0 and B_ij = 0
    [ABD, FMLm, FMLb] = ABD_Matrix_Generator(Al, GF, fml);
    
    % Laminate Modulus
    switch Modulus
        case 'Youngs'
            % Young's Modulus
            E_x     = FMLm.Ex;
            E_y     = FMLm.Ey;
            v_xy    = FMLm.vxy;
            v_yx    = FMLm.vyx;
            G_xy    = FMLm.Gxy;
            E       = E_x;
            v       = v_xy;
        case 'Flexural'
            % Flexural Modulus
            E_x     = FMLb.E1;
            E_y     = FMLb.E2;
            v_xy    = FMLb.v12;
            v_yx    = FMLb.v21;
            G_xy    = FMLb.G12;
            E       = E_x;
            v       = v_xy;
    end
end

% Membrane stiffness
AExx1 = E*t;
AExx0 = E*t*2 + E_a*t_a;

% Bending stifness
EIxx1 = E*t^3/12;
EIxx0 = 2*E*(t^3/12+t*(t/2+t_a/2)^2)+E_a*t_a^3/12;

% Element length
dla = l_A0/q;
dlb = l_B0/q;

%% Module 3: Discretize adherents

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

%% Module 4: Overlap Edge Loads

% Overlap edge loads (minumum and maximum)
[M, Q] = Overlap_Edge_Loads(x1, x0, P, EIxx1, EIxx0, l_A, l_B, t, t_a, BC);

% Extract solutions
M_k  = M.A; % Overlap bending moment based on M_1(x1=l_A)
Q_k  = Q.A; % Overlap shear force based on Q_1(x1=l_A)
Q_kc = Q.C; % Overlap shear force based on equilibrium
M_0 = M.B;  % Overlap bending moment based on M_0(x0=0)
Q_0 = Q.B;  % Overlap shear force based on Q_0(x0=0)

%% Module 5: Adhesive Stresses and Overlap Adherent Load Distributions

% x0 vector must be adjusted to the x-axis system used in the adhesive
% stress analysis: -l_B <= x0 <= 0
x00 = x0-l_B;
x00 = x00.*triu(ones(size(x00)));

% Adhesive stresses derrived using the horizontal force component
F = P*cos(alpha);

% Luo and Tong (2004, 2007) > adhesive thickness inlcuded
[Shear_a, Peel_a, Loads_ad1, Loads_ad2] = Adhesive_Stresses_Adherent_Loads(x00, F, M_k, Q_kc, EIxx1, AExx1, l_B, E, t, E_a, G_a, t_a);

% Goland and Reissner (1944) > adhesive thickness excluded (ta = 0)
%[Shear_a2, Peel_a2] = Adhesive_Stresses_GR(x00, P, l_B, E, t, Ea, Ga, ta, v);

%% Module 6: AL Facesheet Strain and Stress Cycle

if exist('Al','var') && exist ('GF','var') && exist ('fml','var')
    % Stress and strain cycles in the AL facesheets (adjecent to the
    % adhesive interface)
    [Sxx, exx] = Stress_Strain_Cycle(Loads_ad1.N, Loads_ad1.M, Loads_ad2.N, Loads_ad2.M, ABD, AExx1, EIxx1);
    
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
end

%% Module 7: Strain Energy Release Rate

% SERR and Mode Ratio components - Specimen with finite length
[serr, mr] = Strain_Energy_Release_Rate(F, M_k, M_0, E, t, AExx1, EIxx1, 'Fern1und1991');

% SERR and Mode Ratio componentes - Specimen with infinite length
[serr2, mr2] = Strain_Energy_Release_Rate(F, M_k, M_0, E, t, AExx1, EIxx1, 'Brussat');

% Extract solutions
G   = serr.G;
GI  = serr.GI;
GII = serr.GII;
MR  = mr;

% Enforce DAF effect
q1      = round(1/3*q_max);
q2      = round(2/3*q_max);
daf_I   = 0.75;
daf_II  = 0;
x_daf   = 0:1:(q2-q1);

DAF_I   = daf_I*sin(pi/x_daf(end)*x_daf)';
DAF_II  = daf_II*sin(pi/x_daf(end)*x_daf)';

GI(q1:q2,:,:) = GI(q1:q2,:,:)-DAF_I.*GI(q1:q2,:,:);
GII(q1:q2,:,:) = GII(q1:q2,:,:)-DAF_II.*GII(q1:q2,:,:);
G = GI + GII;

MR = GII./G;

%% Module 8: Disbond growth rate and Load Cycles

% Crack Growth Rate (model from D. Burger (2005), FM94 adhesive)
[dbdN, dG1_eq] = Crack_Growth_Rate(GI, GII, MR, c_0, m_0, c_100);

% Load cycles for the disbond to grow the element distance
nraw = dlb./dbdN;
% Round and store removed decimals
n = floor(nraw);
ncut = nraw - n;

%% Module 9: Fatigue Accumulation

% Transform amplitude to a Sm = 0 cycle using the Goodman Relation
Sa_nom_sm0_ad1 = Sa_nom_ad1./(1+Sm_nom_ad1/Al.Su);
Sa_nom_sm0_ad2 = Sa_nom_ad2./(1+Sm_nom_ad2/Al.Su);

% Temp value for now
R_SN = 0.1;

% Find equivalent S-N amplitude stress cycle
Sa_SN_ad1 = Sa_nom_sm0_ad1./(1+Sa_nom_sm0_ad1*((2/(1-R_SN)-1)/Al.Su));
Sa_SN_ad2 = Sa_nom_sm0_ad2./(1+Sa_nom_sm0_ad2*((2/(1-R_SN)-1)/Al.Su));

% Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
% Page 115; Aluminium 2024; Bare sheet, 0.090-inch (2.286 mm) thick
S_eq_ad1 = Sxx.ad1(:,:,2).*(1-R_nom_ad1).^0.56; % Equivalent amplitude S
S_eq_ad1 = S_eq_ad1*1.45038e-7;                 % Pa to ksi
a = S_eq_ad1-15.8;
a(a<0) = 0;
Nf_ad1 = 10.^(11.1-3.97*log10(a));              % Nr. cycles untill fatigue initiation

% Military Handbook - Metallic Materials and Elements for Aerospace Vehicle Structures
% Page 111; Aluminium 2024; Rolled bar
S_eq_ad2 = Sa_nom_ad1.*(1-R_nom_ad1).^0.52;
S_eq_ad2 = S_eq_ad2*1.45038e-7;
Nf_ad2 = 10.^(20.83-9.09*log10(S_eq_ad2));

% Fatigue damage accumulation
% > For now use the rolled bar; the -15.8 in the thin sheet causes issues;
% <-15.8 is the fatigue limit?
Minor = n./Nf_ad1;                              % Damage per cracked element
MinorSum = cumsum(Minor,1);                     % Accumulated total damage

%% Module 10: Results Plotting

% Store in application data to allow acces by the GUI
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
setappdata(0,'GI', GI);
setappdata(0,'GII', GII);
setappdata(0,'dG1_eq', dG1_eq);
setappdata(0,'x',x00*1000);
setappdata(0,'Sy',Al.y1);

run('GUI_FinalPlots');
run('GUI_FinalPlots2');
