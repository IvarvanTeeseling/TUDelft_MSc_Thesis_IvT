%% Matlab reset

clear;
close all;
clc;

%% Module 0: Input

tic

% Adherent material
AdherentSelect = {'Johnson 1986' 'GLARE 4/3 0.4' 'GLARE 4/3 0.3' 'GLARE 5/4 0.4' 'CLT Verification'};
AdherentSelect = AdherentSelect{2}
switch AdherentSelect
    
    case 'Johnson 1986'
        % Data from Two-dimensional analytical solution of elastic stresses
        % for balanced single-lap joints�Variational method (2014)
        % Adherent type indicator
        AdMat   = 'Metal';
        
        % Monolithic adherent properties
        E   = 72.45e9;
        t   = 0.003175; % [m]
        v   = 0.33; % [-]
    case 'GLARE 4/3 0.4'
        % Adherent type indicator
        AdMat   = 'FML';
        
        % Aluminum ply properties
        Al.E1   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.E2   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.G    = 27068e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.t    = 0.4e-3;           % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.Su   = 469e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y1   = 347e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y2   = 299e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.ct1  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.ct2  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Glass Fibre ply properties
        GF.E1   = 48900e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.E2   = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.G    = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.t    = 0.133e-3;         % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct1  = 6.1e-6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct2  = 26.2e-6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Lay-up (1 = AL, 2 = GF)
        fml.layer   = [1 2 2 1 2 2 1 2 2 1];
        fml.theta   = [0 0 0 0 0 0 0 0 0 0];
        fml.t       = (fml.layer == 1)*Al.t + (fml.layer == 2)*GF.t;
        t = sum(fml.t);
        
    case 'GLARE 4/3 0.3'
        % Adherent type indicator
        AdMat   = 'FML';
        
        % Aluminum ply properties
        Al.E1   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.E2   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.G    = 27068e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.t    = 0.3e-3;           % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.Su   = 469e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y1   = 347e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y2   = 299e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.ct1  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.ct2  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Glass Fibre ply properties
        GF.E1   = 48900e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.E2   = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.G    = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.t    = 0.133e-3;         % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct1  = 6.1e-6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct2  = 26.2e-6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Lay-up (1 = AL, 2 = GF)
        fml.layer   = [1 2 2 1 2 2 1 2 2 1];
        fml.theta   = [0 0 0 0 0 0 0 0 0 0];
        fml.t       = (fml.layer == 1)*Al.t + (fml.layer == 2)*GF.t;
        t = sum(fml.t);
        
    case 'GLARE 5/4 0.4'
        % Adherent type indicator
        AdMat   = 'FML';
        
        % Aluminum ply properties
        Al.E1   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.E2   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.G    = 27068e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.t    = 0.4e-3;           % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.Su   = 469e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y1   = 347e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y2   = 299e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.ct1  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.ct2  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Glass Fibre ply properties
        GF.E1   = 48900e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.E2   = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.G    = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.t    = 0.133e-3;         % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct1  = 6.1e-6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct2  = 26.2e-6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Lay-up (1 = AL, 2 = GF)
        fml.layer   = [1 2 2 1 2 2 1 2 2 1 2 2 1];
        fml.theta   = [0 0 0 0 0 0 0 0 0 0 0 0 0];
        fml.t       = (fml.layer == 1)*Al.t + (fml.layer == 2)*GF.t;
        t = sum(fml.t);
        
    case 'CLT Verification'
        % Adherent type indicator
        AdMat   = 'FML';
        
        % Aluminum ply properties
        Al.E1   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.E2   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.G    = 27068e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.t    = 0.4e-3;           % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.Su   = 469e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y1   = 347e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.y2   = 299e6;            % R.C. Alderliesten. Fatigue crack growth propagation and delamination growth in GLARE. PhD thesis, TU Delft, May 2005
        Al.ct1  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.ct2  = 2.32e-5;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Glass Fibre ply properties
        GF.E1   = 48900e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.E2   = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.G    = 5500e6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.t    = 0.133e-3;         % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct1  = 6.1e-6;           % Data fom: Laminate Stiffness Calculator version rca.xls
        GF.ct2  = 26.2e-6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        
        % Lay-up (1 = AL, 2 = GF)
        fml.layer   = [1 2 1 2 1 2 1];
        fml.theta   = [0 0 0 90 0 0 0];
        fml.t       = (fml.layer == 1)*Al.t + (fml.layer == 2)*GF.t;
        t = sum(fml.t);
end

% Adhesive material
AdhesiveSelect = {'FM94 K.06' 'GR_Verification' 'Adh Load Verification' 'Johnson 1986'};
AdhesiveSelect = AdhesiveSelect{1}
switch AdhesiveSelect
    case 'FM94 K.06'
        t_a      = 0.133*10^(-3);                   % Standard thickness
        G_a      = 394e6;                           % Data fom: Laminate Stiffness Calculator version rca.xls
        v_a      = 0.27;                            % Data fom: Laminate Stiffness Calculator version rca.xls
        E_a      = G_a*(2*(1+v_a));                 % Assumed
        c_0      = 5.27*10^(-17);                   % D. Burger PhD report, page 71
        m_0      = 3.78;                            % D. Burger PhD report, page 71
        c_100    = 10^(-17.68379446640316);         % D. Burger PhD report, page 72, figure 5.11, data read using: https://apps.automeris.io/wpd/
    case 'GR_Verification'
        % Data from Modeling of Adhesively Bonded Joints, page 43
        t_a  = 0.078*0.0016;
        E_a  = 70e9*0.04;
        G_a  = 1.1e9;
        v_a  = 0.4;
        % TEMP > not real data
        c_0      = 5.27*10^(-17);
        m_0      = 3.78;
        c_100    = 10^(-17.68379446640316);
    case 'Adh Load Verification'
        % Data from Two-dimensional analytical solution of elastic stresses
        % for balanced single-lap joints�Variational method (2014)
        t_a  = 0.2e-3;  % [m]
        E_a  = 0.7e9; % [Pa]
        G_a  = 0.7e9; % [Pa]
        v_a  = 0.4;
        % TEMP > not real data
        c_0      = 5.27*10^(-17);
        m_0      = 3.78;
        c_100    = 10^(-17.68379446640316);
    case 'Johnson 1986'
        % CLS-A from ASTM Round Robin by Johnson 1986
        t_a  = 0.13e-3;  % [m]
        E_a  = 1.932e9; % [Pa]
        G_a  = E_a/(2*(1+0.4)); % [Pa]
        v_a  = 0.4;
        % TEMP > not real data
        c_0      = 5.27*10^(-17);
        m_0      = 3.78;
        c_100    = 10^(-17.68379446640316);
end

% CLS configuration
ConfigSelect = {'BOPACS' 'WSLS R. Hanx' 'Run_00_Trials' 'GR_Verification' 'Adh Load Verification' 'Johnson 1986' '01' '02'};
ConfigSelect = ConfigSelect{7}
switch ConfigSelect
    case 'BOPACS'
        % Dimensions from Bachelorthesis K. Hoidus, page 35
        L_1 = 0.05+0.145;
        L_2 = 0.145;
        d   = 0.05;
        b_0 = 0;
    case 'WSLS R. Hanx'
        % Dimensions of the WSLS from R. Hanx thesis
        L_1 = 0.400-0.030;
        L_2 = 0.030;
        d   = 0.500;
        b_0  = 0;
    case 'Run_00_Trials'
        % Run_00_Trials CLS configuration
        L_1 = 0.07+0.145;
        L_2 = 0.145;
        d   = 0.048;
        b_0 = 0.0254;
    case 'GR_Verification'
        % Data from Modeling of Adhesively Bonded Joints, page 43
        L_1 = 0.01+0.05;
        L_2 = 0.01;
        d   = 0.018;
        b_0 = 0;
    case 'Adh Load Verification'
        L_1 = 0.04+0.01;
        L_2 = 0.01;
        d   = 0.05;
        b_0 = 0;
    case 'Johnson 1986'
        L_1 = 0.305;
        L_2 = 0.254;
        d   = 0.0254;
        b_0 = 0.0;
    case '01'
        % Run_01_CLS CLS configuration
        L_1 = 0.06835+0.145;
        L_2 = 0.145;
        d   = 0.04366;
        b_0 = 0.0;
    case '02'
        % Run_02_CLS CLS configuration
        L_1 = 0.06835+0.110;
        L_2 = 0.110;
        d   = 0.04366;
        b_0 = 0.03;
end

% Plane stress/strain
StressState = {'Plane Stress' 'Plane Strain'};
StressState = StressState{1}

% Applied load
LCSelect = {'BOPACS K. Hoidus' 'WSLS R. Hanx' 'Johnson 1986' 'Run_00_Trials' 'Run_01_CLS_30%' 'Run_01_CLS_30kN' 'Run_01_CLS_26kN'};
LCSelect = LCSelect{7}
switch LCSelect
    case 'BOPACS K. Hoidus'
        % Load case from: 100% from the load in Bachelorthesis K. Hoidus P49
        P_max   = 60*1e3/d;
        R_load  = 0.1;
    case 'WSLS R. Hanx'
        % Load case from R. Hanx
        P_max   = 1.5*230e3;
        R_load  = 0.1;
    case 'Johnson 1986'
        P_max   = 4.378e5;
        R_load  = 0.1;
    case 'Run_00_Trials'
        P_max   = 24.445e3/d;
        R_load  = 0.1;
    case 'Run_01_CLS_30%'
        P_max   = 22.0327e3/d;
        R_load  = 0.1;
    case 'Run_01_CLS_30kN'
        P_max   = 30e3/d;
        R_load  = 0.1;
    case 'Run_01_CLS_26kN'
        P_max   = 26e3/d;
        R_load  = 0.1;
end

% Support boundary conditions
BC = {'RR' 'CC'};
BC = BC{2}

% Numerical settings - Number of elements
q = 2000;
% Numerical settings - Number of cracked elements
q_max = 1500;
% Numerical settings - Discretization method
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{2}

%% Module 1: Load Paramaters

% Applied force
P_min    = P_max*R_load;    % [N/m]
P(1,1,1) = P_min;           % [N/m]
P(1,1,2) = P_max;           % [N/m]

% Applied stress
S_max    = P_max/t;         % [N/m^3]
S_min    = S_max*R_load;    % [N/m^3]
S(1,1,1) = S_min;           % [N/m^3]
S(1,1,2) = S_max;           % [N/m^3]

%% Module 2: FML Laminate Properties

if strcmp(AdMat, 'FML')
    % ABD matrix and equivalent laminate properties
    [ABD, FML] = ABD_Matrix_Generator(Al, GF, fml, 'mem');
    
    % Isolate output
    E_x     = FML.Ex;
    E_x_ps  = FML.Ex_ps;
    E_y     = FML.Ey;
    v_xy    = FML.vxy;
    v_yx    = FML.vyx;
    G_xy    = FML.Gxy;
    E       = E_x;
    v       = v_xy;
end

% Adjust for plane strain if applicable
if strcmp(StressState, 'Plane Strain')
    if strcmp(AdMat, 'FML')
        E = E_x_ps;
    else
        E = E/(1-v^2);
    end
end

% Membrane stiffness (only for identical adherends)
AExx_A       = E*t;
AExx_B       = E*t*2 + E_a*t_a;
AExx_B_ta0   = E*t*2;                                        % Assuming t_a~0

% Bending stifness (only for identical adherends)
EIxx_A       = E*t^3/12;
EIxx_B       = 2*E*(t^3/12+t*(t/2+t_a/2)^2)+E_a*t_a^3/12;
EIxx_B_ta0   = E*(2*t)^3/12;                                 % Assuming t_a~0

%% Module 3: Discretize adherends

% Initital geometry
l_A0 = L_1-(L_2-b_0);
l_B0 = L_2-b_0;

alpha = (t+t_a)./(2*(l_A0+l_B0));

% Discretize adherents according to the xA, xB and xBB reference frame
[xA, xAB, xB, xBB, l_A, l_B] = Discretize_Adherents(l_A0, l_B0, q, q_max, NaN, DiscretizeMethod);

%% Module 4: Overlap Edge Loads

% Overlap edge loads (minumum and maximum)
[M_k, M_k0, Q_k, Q_k0, V_k, M, Q, w] = Overlap_Edge_Loads(xA, xB, P, EIxx_A, EIxx_B_ta0, l_A, l_B, t, 0, BC);

%% Module 5: Adhesive Stresses and Overlap Adherent Load Distributions

% Adhesive stresses derrived using the horizontal force component
F = P*cos(alpha);

% Luo and Tong (2004, 2007) > adhesive thickness inlcuded
[Shear_a, Peel_a] = Overlap_Adhesive_Stresses(xBB, F, M_k, V_k, l_B, E, t, E_a, G_a, t_a);

% Get the N, M and Q distributions in the overlap region
[OALoad_T, OALoad_L] = Overlap_Adherent_Load_Distributions(xB, t, t_a, F, V_k, M_k, Shear_a, Peel_a, 'num');

toc

%% Module 6: AL Facesheet Strain and Stress Cycle

if strcmp(AdMat, 'FML')
    % Free adherent
    [~, e_xx_A, ~, Sm_nom_A, Sa_nom_A] = Stress_Strain_Cycle(P, M.A, ABD, AExx_A, EIxx_A);
    
    % Overlap region - Top adherent
    [~, e_xx_BT, ~, Sm_nom_BT, Sa_nom_BT] = Stress_Strain_Cycle(OALoad_T.N, OALoad_T.M, ABD, AExx_A, EIxx_A);
    
    % Overlap region - Lower adherent
    [~, e_xx_BL, ~, Sm_nom_BL, Sa_nom_BL] = Stress_Strain_Cycle(OALoad_L.N, OALoad_L.M, ABD, AExx_A, EIxx_A);
end

%% Module 7: Strain Energy Release Rate

% Strain Energy Release Rate (3 methods)
[serr] = Strain_Energy_Release_Rate('Fern1und1991', AExx_A, EIxx_A, P, M_k, M_k0);

% Load the SERR effect
%load('SERR_DIFF_1_vs_2.mat');
%load('SERR_DIFF_Fc100.mat');
load('SERR_DIFF_Fc130.mat');

% Account for the DAF effect
[serr] = DAF_SERR_Effect(xB(1,1:q_max+1)', serr.GI, serr.GII, serr_diff_Fc130(5:end,1,2)/1000-0.045, serr_diff_Fc130(5:end,4,2)/100, serr_diff_Fc130(5:end,3,1)/100);
[serr2] = DAF_SERR_Effect(xB(1,1:q_max+1)', serr.GI, serr.GII, serr_diff_Fc130(5:end,1,1)/1000-0.045, serr_diff_Fc130(5:end,4,1)/100, serr_diff_Fc130(5:end,3,1)/100);
%[serr] = DAF_SERR_Effect(xB(1,1:q_max+1)', serr.GI, serr.GII, serr_diff_Fc100(5:end,1,2)/1000-0.045, serr_diff_Fc100(5:end,4,2)/100, zeros(size(serr_diff_Fc100(5:end,3,2))));

%% Module 8: Adhesive DGR and numerical integration

% Crack Growth Rate (model from D. Burger (2005), FM94 adhesive)
[dbdN, dG1_eq] = Crack_Growth_Rate(serr.GI, serr.GII, serr.MR, c_0, m_0, c_100, 1e-30);

if dbdN(1) == 0
    disp('Warning! No adhesive disbond growth exists..!');
end

% Set crack increments to stop if dbdN = 0 has occured
b = xB(1, 1:q_max)';
if any(dbdN==0)
    ind = find(dbdN==0, 1);
    b(ind:end) = b(ind);
    b = b+b_0;
end

% dN   = l_B0/q./dbdN;
% dN1  = l_B0/q./dbdN1;
% dN2  = l_B0/q./dbdN2;
%
% figure(1)
% hold on
% plot(cumsum(dN(1:end-1)), b*1000, 'r')
% plot(cumsum(dN1(1:end-1)), b*1000, 'b')
% plot(cumsum(dN2(1:end-1)), b*1000, 'g')
% plot(cumsum(dN2(1:end-1)), ones(size(b))*50, '-c')
% hold off
% grid on
% xlabel('Load Cycle (N) [-]')
% ylabel('Disbond Length (b) [mm]')
% legend('No DAF', 'Clamped DAF', 'Clamped DAF - Neutral % G_{II} Footprint', 'Clamped DAF position')
% title('Arrested Disbond Growth Behavior (F_t = 26kN, F_c = 10N/mm^2)')
% grid on


%% Module 9: Fatigue Accumulation

if strcmp(AdMat, 'FML')
    % First we have to merge the stress matrices of region A and B of the
    % upper adherent so that we get one stress matrix spanning all elements
    [Sa_nom_ACB] = xMatrix_Zipper(Sa_nom_A, Sa_nom_BT, q, q_max, NaN, 'zip');
    [Sm_nom_ACB] = xMatrix_Zipper(Sm_nom_A, Sm_nom_BT, q, q_max, NaN, 'zip');
    
    % Fatigue damage accumulation
    [Minor, dMinor, dN] = Adherent_Fatigue_Accumulation('Military Handbook - Sheet', Sa_nom_ACB(:,:,:), Sm_nom_ACB(:,:,:), Al.Su, dbdN, l_B0/q);
end

%% Module 10: Results Plotting

figures = 0;

if figures == 1
    
    n1 = 1;
    n2 = 1200;
    
    figure(1)
    hold on
    plot([0 0], [-1e-3 6e-3], 'g')
    plot([xA(n1,:)-l_A(n1) xB(n1,:)]*1000, [e_xx_A(n1,:,2,2) e_xx_BT(n1,:,2,2)],'b')
    plot([xA(n1,:)-l_A(n1) xB(n1,:)]*1000, [e_xx_A(n1,:,2,1) e_xx_BT(n1,:,2,1)],'--b')
    plot(xB(n1,:)*1000, e_xx_BL(n1,:,2,2),'--r')
    plot(xB(n1,:)*1000, e_xx_BL(n1,:,2,1),'r')
    hold off
    grid on
    legend('Adhesive crack','Top FML (top ply)','Top FML (bottom ply)','Bottom FML (top ply)','Bottom FML (bottom ply)')
    xlabel('x-position [mm]')
    ylabel('\epsilon_{xx} [-]')
    
    figure(2)
    hold on
    plot(xA(n1,:)-l_A(n1), M.A(n1,:,2))
    plot(xB(n1,:), M.B(n1,:,2))
    hold off
    
    figure(3)
    hold on
    plot(xA(n1,:), w.A(n1,:,2),'r')
    plot(xB(n1,:)+l_A(n1), w.B(n1,:,2),'b')
    plot(xA(n2,:), w.A(n2,:,2),'--r')
    plot(xB(n2,:)+l_A(n2), w.B(n2,:,2),'--b')
    hold off
    
    figure(4)
    hold on
    plot(xA(n2,:)-xA(n2,end), Q.A(n2,:,2))
    plot(xB(n2,:), Q.B(n2,:,2))
    hold off
    
    x = -24:0.1:24;
    x = repmat(x,length(xB(n2,:)),1);
    y = repmat(xB(n2,:)',1,size(x,2));
    z = repmat(e_xx_BL(n2,:,2,1)',1,size(x,2));
    
    figure(5)
    contourf(x,y,z, linspace(min(z(:)), max(z(:)), 15))
    colorbar
    caxis([min(z(:)) max(z(:))])
    
end

gui = 1;

if gui == 1
    % Store in application data to allow acces by the GUI
    if strcmp(AdMat, 'FML')
        setappdata(0, 'l_A', l_A*1000);
        setappdata(0,'N1',OALoad_T.N);
        setappdata(0,'Q1',OALoad_T.Q);
        setappdata(0,'M1',OALoad_T.M);
        setappdata(0,'N2',OALoad_L.N);
        setappdata(0,'Q2',OALoad_L.Q);
        setappdata(0,'M2',OALoad_L.M);
        setappdata(0,'Shear_a',Shear_a);
        setappdata(0,'Peel_a',Peel_a);
        setappdata(0,'N',cumsum(dN));
        setappdata(0,'b',b*1000);
        setappdata(0,'dbdN',dbdN(1:length(dN)));
        setappdata(0,'dMinor',dMinor);
        setappdata(0,'Minor',Minor);
        setappdata(0,'Sm_nom_ACB',Sm_nom_ACB);
        setappdata(0,'Sa_nom_ACB',Sa_nom_ACB);
        setappdata(0,'Sm_nom_BT',Sm_nom_BT);
        setappdata(0,'Sa_nom_BT',Sa_nom_BT);
        setappdata(0,'Sm_nom_BL',Sm_nom_BL);
        setappdata(0,'Sa_nom_BL',Sa_nom_BL);
        setappdata(0,'Sm_nom_A',Sm_nom_A);
        setappdata(0,'Sa_nom_A',Sa_nom_A);
        setappdata(0,'GI', serr.GI);
        setappdata(0,'GII', serr.GII);
        setappdata(0,'dG1_eq', dG1_eq);
        setappdata(0,'xAB',xAB*1000);
        setappdata(0,'x',xBB*1000);
        setappdata(0,'Sy',Al.y1);
        setappdata(0,'Su',Al.Su);
        
        %run('GUI_FinalPlots');
        run('GUI_FinalPlots3');
    end
end

toc