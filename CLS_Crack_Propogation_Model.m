%% Matlab reset

clear;
close all;
clc;

%% Module 0: Input

% Adherent material
AdherentSelect = {'GR Verification' 'Adh Load Verification' 'D. Burger PhD Validation' 'GLARE 4/3 0.4' 'GLARE 4/3 0.3' 'CLT Verification' 'Johnson 1986'};
AdherentSelect = AdherentSelect{7}
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
AdhesiveSelect = AdhesiveSelect{4}
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
ConfigSelect = ConfigSelect{6}
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
LCSelect = {'BOPACS K. Hoidus' 'WSLS R. Hanx' 'Johnson 1986' 'Run_00_Trials'};
LCSelect = LCSelect{3}
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
end

% Support boundary conditions
BC = {'RR' 'CC'};
BC = BC{2}

% Numerical settings - Number of elements
q = 3000;
% Numerical settings - Number of cracked elements
q_max = ceil(q*9/10);
% Numerical settings - Discretization method
DiscretizeMethod = {'LeftBoundary' 'Central' 'RightBoundary'};
DiscretizeMethod = DiscretizeMethod{2}

%% Module 1: Load Paramaters

% Applied force
P_min    = P_max*R_load;    % [N/m]
P(1,1,1) = P_min;           % [N/m]
P(1,1,2) = P_max;           % [N/m]

% Applied stress
S_max    = P_max*t;         % [N/m^2]
S_min    = S_max*R_load;    % [N/m^2]
S(1,1,1) = S_min;           % [N/m^2]
S(1,1,2) = S_max;           % [N/m^2]

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

% Angle neutral axis w.r.t. adherent neutral axis
alpha = (t+t_a)./(2*(l_A0+l_B0));

% Element length
dlA = l_A0/q;
dlB = l_B0/q;

switch DiscretizeMethod
    case 'LeftBoundary'
        xA = 0:dlA:l_A0-dlA;
        xB = 0:dlB:l_B0-dlB;
    case 'Central'
        xA = dlA/2:dlA:l_A0-dlA/2;
        xB = dlB/2:dlB:l_B0-dlB/2;
    case 'RightBoundary'
        xA = dlA:dlA:l_A0;
        xB = dlB:dlB:l_B0;
end

% x1 matrix where x1 spans l_A (free adherent)
xA = xA.*ones(q_max,1);
xA = xA-dlA*ones(q_max,q).*(0:q_max-1)';
xA = xA.*triu(ones(size(xA)));

% x0 matrix where x0 spans l_B (overlap region)
xB = xB.*ones(q_max,1);
xB = xB-dlB*ones(q_max,q).*(0:q_max-1)';
xB = xB.*triu(ones(size(xB)));

% l_A and l_B for each crack increment
l_A = l_A0*ones(q_max,1)+dlB*(0:q_max-1)';
l_B = l_B0*ones(q_max,1)-dlB*(0:q_max-1)';

%% Module 4: Overlap Edge Loads

% Overlap edge loads (minumum and maximum)
[M_k, M_k0, Q_k, Q_k0, V_k, M, Q, w] = Overlap_Edge_Loads(xA, xB, P, EIxx_A, EIxx_B_ta0, l_A, l_B, t, 0, BC);


%% Module 5: Adhesive Stresses and Overlap Adherent Load Distributions

% x0 vector must be adjusted to the x-axis system used in the adhesive
% stress analysis: -l_B <= xB <= 0
xBB = xB-l_B;
xBB = xBB.*triu(ones(size(xBB)));

% Adhesive stresses derrived using the horizontal force component
F = P*cos(alpha);

% Luo and Tong (2004, 2007) > adhesive thickness inlcuded
[Shear_a, Peel_a] = Overlap_Adhesive_Stresses(xBB, F, M_k, V_k, l_B, E, t, E_a, G_a, t_a);

[Loads_ad1, Loads_ad2] = Overlap_Adherent_Load_Distributions(xB, t, t_a, F, V_k, M_k, Shear_a, Peel_a, 'num');

%% Module 6: AL Facesheet Strain and Stress Cycle

if strcmp(AdMat, 'FML')
    % Free adherent
    [~, e_xx_A, ~, ~, ~] = Stress_Strain_Cycle(P, M.A, ABD, AExx_A, EIxx_A);
    
    % Overlap region - Top adherent
    [S_xx_B1, e_xx_B1, R_nom_B1, Sm_nom_B1, Sa_nom_B1] = Stress_Strain_Cycle(Loads_ad1.N, Loads_ad1.M, ABD, AExx_A, EIxx_A);
    
    % Overlap region - Bottom adherent
    [~, e_xx_B2, ~, ~, ~] = Stress_Strain_Cycle(Loads_ad2.N, Loads_ad2.M, ABD, AExx_A, EIxx_A);
end

%% Module 7: Strain Energy Release Rate

% Strain Energy Release Rate (3 methods)
[serr, mr]          = Strain_Energy_Release_Rate('Fern1und1991', AExx_A, EIxx_A, P, M_k, M_k0);
[serr2, mr2]        = Strain_Energy_Release_Rate('Verreman1992', t_a, E_a, G_a, max(Peel_a,[],2), max(Shear_a,[],2));
[serr_inf, mr_inf]  = Strain_Energy_Release_Rate('Brussat1977', t, E, AExx_A, EIxx_A, EIxx_B_ta0, P, M_k, M_k0);

% Enforce DAF effect
q1      = round(1/3*q_max);
q2      = round(2/3*q_max);
daf_I   = 0.8;
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
[dbdN, dG1_eq] = Crack_Growth_Rate(serr.GI, serr.GII, serr.MR, c_0, m_0, c_100, 10^(-8));

if dbdN(1) == 0
    disp('Warning! No adhesive disbond growth exists..!');
end

%% Module 9: Fatigue Accumulation

if strcmp(AdMat, 'FML')
    % Fatigue damage accumulation
    [Minor_csm, Minor, dN, N_f] = Adherent_Fatigue_Accumulation('Military Handbook - Sheet', Sa_nom_B1(:,:,1), Sm_nom_B1(:,:,1), R_nom_B1(:,:,1), Al.Su, dbdN, dlB);
end

%% Module 10: Results Plotting

figures = 1;

if figures == 1
    figure(1)
    hold on
    plot([0 0], [-1e-3 6e-3], 'g')
    plot([xA(1,:)-xA(1,end) xB(1,:)]*1000, [e_xx_A(1,:,2,2) e_xx_B1(1,:,2,2)],'b')
    plot([xA(1,:)-xA(1,end) xB(1,:)]*1000, [e_xx_A(1,:,2,1) e_xx_B1(1,:,2,1)],'--b')
    plot(xB(1,:)*1000, e_xx_B2(1,:,2,2),'--r')
    plot(xB(1,:)*1000, e_xx_B2(1,:,2,1),'r')
    hold off
    grid on
    legend('Adhesive crack','Top FML (top ply)','Top FML (bottom ply)','Bottom FML (top ply)','Bottom FML (bottom ply)')
    xlabel('x-position [mm]')
    ylabel('\epsilon_{xx} [-]')
    
    figure(2)
    hold on
    plot(xA(1,:)-xA(1,end), M.A(1,:,2))
    plot(xB(1,:), M.B(1,:,2))
    hold off
    
    figure(3)
    hold on
    plot(xA(1,:)-xA(1,end), w.A(1,:,2))
    plot(xB(1,:), w.B(1,:,2))
    hold off
    
    figure(4)
    hold on
    plot(xA(1,:)-xA(1,end), Q.A(1,:,2))
    plot(xB(1,:), Q.B(1,:,2))
    hold off
    
    x = -24:0.1:24;
    x = repmat(x,length(xB(1,:)),1);
    y = repmat(xB(1,:)',1,size(x,2));
    z = repmat(e_xx_B2(1,:,2,1)',1,size(x,2));
    
    figure(3)
    contourf(x,y,z, linspace(min(z(:)), max(z(:)), 15))
    colorbar
    caxis([min(z(:)) max(z(:))])
end

gui = 0;

if gui == 1
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
        setappdata(0,'Sm_nom_ad1',Sm_nom_B1(:,:,1));
        setappdata(0,'Sa_nom_ad1',Sa_nom_B1(:,:,1));
        setappdata(0,'GI', serr.GI);
        setappdata(0,'GII', serr.GII);
        setappdata(0,'dG1_eq', dG1_eq);
        setappdata(0,'G_inf', serr_inf.G);
        setappdata(0,'x',xBB*1000);
        setappdata(0,'Sy',Al.y1);
        setappdata(0,'Su',Al.Su);
        
        run('GUI_FinalPlots');
        run('GUI_FinalPlots2');
    end
end