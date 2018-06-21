%% Matlab reset

clear;
close all;
clc;

%% Module 0: Input

% Adherent material
AdherentSelect = {'GR Verification' 'Adh Load Verification' 'D. Burger PhD Validation' 'GLARE 4/3 0.4' 'GLARE 4/3 0.3' 'CLT Verification' 'Johnson 1986'};
AdherentSelect = AdherentSelect{4}
switch AdherentSelect
    case 'GR Verification'
        % Data from Modeling of Adhesively Bonded Joints, page 43
        % Adherent type indicator
        AdMat   = 'Metal';
        
        % Input from E. Verreman
        E   = 70e9;
        t   = 0.0016;
        v   = 0.34;
    case 'Adh Load Verification'
        % Data from Two-dimensional analytical solution of elastic stresses
        % for balanced single-lap joints—Variational method (2014)
        % Adherent type indicator
        AdMat   = 'Metal';
        
        % Monolithic adherent properties
        E   = 70e9;
        t   = 0.004; % [m]
        v   = 0.34; % [-]
    case 'D. Burger PhD Validation'
        % Data from the PhD by Daniel Burger, Chapter 6
        AdMat   = 'Metal';
        
        % Monolithic adherent properties
        E   = 72.45e9;
        t   = 0.006; % [m]
        v   = 0.33; % [-]
    case 'Johnson 1986'
        % Data from Two-dimensional analytical solution of elastic stresses
        % for balanced single-lap joints—Variational method (2014)
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
        Al.Su   = 469e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
        Al.y1   = 347e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
        Al.y2   = 299e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
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
        Al.t    = 0.4e-3;           % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.Su   = 469e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
        Al.y1   = 347e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
        Al.y2   = 299e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
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
        
    case 'CLT Verification'
        % Adherent type indicator
        AdMat   = 'FML';
        
        % Aluminum ply properties
        Al.E1   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.E2   = 72000e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.G    = 27068e6;          % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.v12  = 0.33;             % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.t    = 0.4e-3;           % Data fom: Laminate Stiffness Calculator version rca.xls
        Al.Su   = 469e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
        Al.y1   = 347e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
        Al.y2   = 299e6;            % http://asm.matweb.com/search/SpecificMaterial.asp?bassnum=ma2024t4
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
        E_a      = 1e9;                             % Data fom: Laminate Stiffness Calculator version rca.xls
        G_a      = 394e6;                           % Data fom: Laminate Stiffness Calculator version rca.xls
        v_a      = 0.27;                            % Data fom: Laminate Stiffness Calculator version rca.xls
        
        E_a      = G_a*(2*(1+v_a));
        
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
        % for balanced single-lap joints—Variational method (2014)
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
ConfigSelect = {'BOPACS' 'WSLS R. Hanx' 'Run_00_Trials' 'GR_Verification' 'Adh Load Verification' 'Johnson 1986' 'D. Burger PhD Validation'};
ConfigSelect = ConfigSelect{3}
switch ConfigSelect
    case 'BOPACS'
        % Dimensions from Bachelorthesis K. Hoidus, page 35
        L_1 = 0.05+0.145;
        L_2 = 0.145;
        d   = 0.05;
        b_0  = 0;
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
        d   = 0.05;
        b_0  = 0;
    case 'GR_Verification'
        % Data from Modeling of Adhesively Bonded Joints, page 43
        L_1 = 0.01+0.05;
        L_2 = 0.01;
        d   = 0.05;
        b_0  = 0;
    case 'Adh Load Verification'
        L_1 = 0.04+0.01;
        L_2 = 0.01;
        d   = 0.05;
        b_0  = 0;
    case 'Johnson 1986'
        L_1 = 0.305;
        L_2 = 0.254;
        d   = 0.0254;
        b_0 = 0.101;
    case 'D. Burger PhD Validation'
        % Data from the PhD by Daniel Burger, Chapter 6
        L_1 = 0.300-0.025;
        L_2 = 0.075;
        d   = 0.025;
        b_0 = 0.015;
end

% Plane stress/strain
StressState = {'Plane Stress' 'Plane Strain'};
StressState = StressState{2}

% Applied load
LCSelect = {'BOPACS K. Hoidus' 'WSLS R. Hanx' 'Johnson 1986' 'D. Burger PhD Validation'};
LCSelect = LCSelect{2}
switch LCSelect
    case 'BOPACS K. Hoidus'
        % Load case from: 200% from the load in Bachelorthesis K. Hoidus P49
        S_max   = 1*23.16*1e3/(d*t); % [N/m^2]
        R_load  = 0.1;
    case 'WSLS R. Hanx'
        % Load case from R. Hanx
        S_max   = 1.5*230e3/t;
        R_load  = 0.1;
    case 'Johnson 1986'
        S_max   = 4.378e5/t;
        R_load  = 0.1;
    case 'D. Burger PhD Validation'
        % Data from the PhD by Daniel Burger, Chapter 6
        S_max   = 14e3/(d*t);
        R_load  = 0.1;
end

% Support boundary conditions (only for overlap loads)
BC = {'RR' 'CC'};
BC = BC{2}

% Numerical settings - Number of elements
q = 3000;
% Numerical settings - Number of cracked elements
q_max = ceil(q*9/10);
% Numerical settings - Discretization method
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{2}

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

if strcmp(AdMat, 'FML')
    
    % ABD matrix and equivalent laminate properties
    [ABD, FML] = ABD_Matrix_Generator(Al, GF, fml, 'mem');
    
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
AExx1       = E*t;
AExx0       = E*t*2 + E_a*t_a;
AExx0_ta0   = E*t*2;                                        % Assuming t_a~0

% Bending stifness (only for identical adherends)
EIxx1       = E*t^3/12;
EIxx0       = 2*E*(t^3/12+t*(t/2+t_a/2)^2)+E_a*t_a^3/12;
EIxx0_ta0   = E*(2*t)^3/12;                                 % Assuming t_a~0

% Element length
dla = l_A0/q;
dlb = l_B0/q;

%% Module 3: Discretize adherends

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

% x1 matrix where x1 spans l_A (free adherent)
x1 = x1.*ones(q_max,1);
x1 = x1-dla*ones(q_max,q).*(0:q_max-1)';
x1 = x1.*triu(ones(size(x1)));

% x0 matrix where x0 spans l_B (overlap region)
x0 = x0.*ones(q_max,1);
x0 = x0-dlb*ones(q_max,q).*(0:q_max-1)';
x0 = x0.*triu(ones(size(x0)));

% l_A and l_B for each crack increment
l_A = l_A0*ones(q_max,1)+dlb*(0:q_max-1)';
l_B = l_B0*ones(q_max,1)-dlb*(0:q_max-1)';

%% Module 4: Overlap Edge Loads

% Overlap edge loads (minumum and maximum)
[M, Q, V] = Overlap_Edge_Loads(x1, x0, P, EIxx1, EIxx0_ta0, l_A, l_B, t, 0, BC);

% Extract solutions
M_k  = M.A; % Overlap bending moment based on M_1(x1=l_A)
Q_k  = Q.A; % Overlap shear force based on Q_1(x1=l_A)
V_k  = V;   % Overlap shear force based on equilibrium
M_0  = M.B; % Overlap bending moment based on M_0(x0=0)
Q_0  = Q.B; % Overlap shear force based on Q_0(x0=0)

%% Module 5: Adhesive Stresses and Overlap Adherent Load Distributions

% x0 vector must be adjusted to the x-axis system used in the adhesive
% stress analysis: -l_B <= x0 <= 0
x00 = x0-l_B;
x00 = x00.*triu(ones(size(x00)));

% Adhesive stresses derrived using the horizontal force component
F = P*cos(alpha);

% Luo and Tong (2004, 2007) > adhesive thickness inlcuded
[Shear_a, Peel_a] = Overlap_Adhesive_Stresses(x00, F, M_k, V_k, l_B, E, t, E_a, G_a, t_a);

% Goland and Reissner (1944) > adhesive thickness excluded (ta = 0)
%[Shear_a2, Peel_a2] = Adhesive_Stresses_GR(x00, P, l_B, E, t, E_a, G_a, t_a, v);
%[Shear_a3, Peel_a3] = Overlap_Adhesive_Stresses_LT(x0, F, M_k, V_k, l_B, E, G, t, E_a, G_a, t_a, AExx1, EIxx1);

[Loads_ad1, Loads_ad2] = Overlap_Adherent_Load_Distributions(x0, t, t_a, F, V_k, M_k, Shear_a, Peel_a, 'num');

%% Module 6: AL Facesheet Strain and Stress Cycle

if strcmp(AdMat, 'FML')
    
    [S_xx, e_xx, R_nom, Sm_nom, Sa_nom] = Stress_Strain_Cycle(Loads_ad1.N, Loads_ad1.M, Loads_ad2.N, Loads_ad2.M, ABD, AExx1, EIxx1);

end

%% Module 7: Strain Energy Release Rate

[serr, mr]          = Strain_Energy_Release_Rate('Fern1und1991', AExx1, EIxx1, P, M_k, M_0);
[serr2, mr2]        = Strain_Energy_Release_Rate('Verreman1992', t_a, E_a, G_a, max(Peel_a,[],2), max(Shear_a,[],2));
[serr_inf, mr_inf]  = Strain_Energy_Release_Rate('Brussat1977', t, E, AExx1, EIxx1, EIxx0_ta0, P, M_k, M_0);

% Enforce DAF effect
q1      = round(1/3*q_max);
q2      = round(2/3*q_max);
daf_I   = 0;
daf_II  = 0;
x_daf   = 0:1:(q2-q1);

DAF_I   = daf_I*sin(pi/x_daf(end)*x_daf)';
DAF_II  = daf_II*sin(pi/x_daf(end)*x_daf)';

serr.GI(q1:q2,:,:)   = serr.GI(q1:q2,:,:)-DAF_I.*serr.GI(q1:q2,:,:);
serr.GII(q1:q2,:,:)  = serr.GII(q1:q2,:,:)-DAF_II.*serr.GII(q1:q2,:,:);
serr.G               = serr.GI + serr.GII;
serr.MR              = serr.GII./serr.G;

%% Module 8: Adhesive DGR and numerical integration

% Crack Growth Rate (model from D. Burger (2005), FM94 adhesive)
[dbdN, dG1_eq] = Crack_Growth_Rate(serr.GI, serr.GII, serr.MR, c_0, m_0, c_100, 10^(-10));

if dbdN(1) == 0
    error('Error. No adhesive disbond growth exists..!');
end

%% Module 9: Fatigue Accumulation

if strcmp(AdMat, 'FML')
        
    [Minor_csm, Minor, dN, N_f] = Adherent_Fatigue_Accumulation('Military Handbook - Sheet', Sa_nom.ad1, Sm_nom.ad1, R_nom.ad1, Al.Su, dbdN, dlb);
    
end

%% Module 10: Results Plotting

% Store in application data to allow acces by the GUI
if strcmp(AdMat, 'FML')
    setappdata(0,'N1',Loads_ad1.N);
    setappdata(0,'Q1',Loads_ad1.Q);
    setappdata(0,'M1',Loads_ad1.M);
    setappdata(0,'N2',Loads_ad2.N);
    setappdata(0,'Q2',Loads_ad2.Q);
    setappdata(0,'M2',Loads_ad2.M);
    setappdata(0,'Shear_a',Shear_a);
    setappdata(0,'Peel_a',Peel_a);
    setappdata(0,'N',cumsum(dN));
    setappdata(0,'dbdN',dbdN(1:length(dN)));
    setappdata(0,'MinorSum',Minor_csm);
    setappdata(0,'Sm_nom_ad1',Sm_nom.ad1);
    setappdata(0,'Sa_nom_ad1',Sa_nom.ad1);
    setappdata(0,'GI', serr.GI);
    setappdata(0,'GII', serr.GII);
    setappdata(0,'dG1_eq', dG1_eq);
    setappdata(0,'G_inf', serr_inf.G);
    setappdata(0,'x',x00*1000);
    setappdata(0,'Sy',Al.y1);
    
    run('GUI_FinalPlots');
    run('GUI_FinalPlots2');
end