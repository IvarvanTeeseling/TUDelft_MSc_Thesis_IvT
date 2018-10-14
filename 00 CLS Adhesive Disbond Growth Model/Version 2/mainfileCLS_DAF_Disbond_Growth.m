%% Matlab reset

clear;
close all;
clc;
format short;
tic
%% User Input - Define Input Here

% 1) Aluminum material selection
% @matAlSelect
%   > Input options:
%       1 = Aluminium 2024-T3 rolled
matAlSelect = 1;

% 2) Fiber ply material selection
% @matGFSelect
%   > Input options:
%       1 = UD glass fiber FM94 S2 prepreg
%       2 = ADD HERE
matGFSelect = 1;

% 3) Adhesive material selection
% @matAdSelect
%   > Input options:
%       1 = FM94U Film Adhesive (without a carrier)
%       2 = ADD HERE
matAdSelect = 1;

% 4) Laminate layup selection
% @layupSelect
%   > Input options:
%       1 = GLARE 2A-4/3-0.3
%       2 = GLARE 2A-4/3-0.4
%       3 = GLARE 2A-5/4-0.4
layupSelect = 2;

% 5) Load cycle selection
% @lcSelect
%   > Input options:
%       1 = BOPACS CLS Specimen, K. Hoidus [N]
%       2 = Thesis  WSLS Specimen, R. Hanx [N/m]
%       3 = An ASTM Round Robin CLS specimen, Johnson (1986) [N]
%       4 = Thesis CLS, I. van Teeseling [N]
lcSelect = 4;

% 6) Support boundary conditions
% @bcSelect
%   > Input options:
%       'CC' = Clamped-Clamped
%       'RR' = Roller-Roller
%       'CR' = Clamped-Roller (CURRENTLY NOT INCLUDED)
bcSelect = 'CC';

% 7) Select the CLS geometry
% @clsSelect
%   > Input options:
%       1 = An ASTM Round Robin (CLS-A) from Johnson 1986
%       2 = CLS Production series 2
%       3 = CLS Production series 4
clsSelect = 2;

% 8) Discretization
% @mesh.qac
%   > Number of elements for the adherents AC section
% @mesh.qcb
%   > Number of elements for the adherents CB section
% @mesh.qcracked
%   > xCB number of elements for disbond propogation
% @mesh.method
%   > Take f(x) at the left, central or right element x-coordinate
%   > Input options:
%       'L' = Element left boundary value
%       'C' = Element central value
%       'R' = Element right boundary value
mesh.qac = 500;
mesh.qcb = 4000;
mesh.qcrack = 500;
mesh.xpos = 'C';

%% Paramater Initialization - Add Options Here

% 1) Metal ply material
if matAlSelect == 1
    % Aluminium 2024-T3 rolled
    alu = Metal('elastic', 72000e6, 72000e6, 27068e6, 0.33, ...
        'sy1', 347e6, 'sy2', 299e6, 'su1', 469e6, ...
        'su2', 469e6, 'ct1', 2.32e-5, 'ct2', 2.32e-5, ...
        'name', 'aluminium 2024-T3 rolled');
end
disp('>')
disp(['> Metal ply material created: ' alu.name])

% 2) Fiber ply material
if matGFSelect == 1
    % UD glass fiber FM94 S2 prepreg
    gf = FRP('elastic', 48900e6, 5500e6, 5500e6, 0.33, ...
        'ct1', 6.1e-6, 'ct2', 26.2e-6, ...
        'name', 'UD glass fiber S2 FM94 prepreg');
end
disp('>')
disp(['> Fiber ply material created: ' gf.name])

% 3) Adhesive material
if matAdSelect == 1
    % FM94U adhesive film (without a carrier)
    adh = Adhesive('elastic', 394e6*(2*(1+0.27)), 394e6*(2*(1+0.27)), ...
        394e6, 0.33, 5.27e-17, 10^(-17.68379446640316),  3.78, ...
        'name', 'FM94U Film Adhesive');
end
disp('>')
disp(['> Adhesive material created: ' adh.name])

% 4) Fiber Metal Laminate material
if layupSelect == 1
    % GLARE 2A-4/3-0.3
    idx = [1 2 2 1 2 2 1 2 2 1];
    tply  = (idx == 1)*0.3e-3 + (idx == 2)*0.133e-3;
    theta = [0 0 0 0 0 0 0 0 0 0];
    layup = [idx ; tply ; theta];
    disp('>')
    disp('> Laminate selected: GLARE 2A-4/3-0.3')
elseif layupSelect == 2
    % GLARE 2A-4/3-0.4
    idx = [1 2 2 1 2 2 1 2 2 1];
    tply  = (idx == 1)*0.4e-3 + (idx == 2)*0.133e-3;
    theta = [0 0 0 0 0 0 0 0 0 0];
    layup = [idx ; tply ; theta];
    disp('>')
    disp('> Laminate selected: GLARE 2A-4/3-0.4')
elseif layupSelect == 3
    % GLARE 2A-5/4-0.4
    idx = [1 2 2 1 2 2 1 2 2 1 2 2 1];
    tply  = (idx == 1)*0.4e-3 + (idx == 2)*0.133e-3;
    theta = [0 0 0 0 0 0 0 0 0 0 0 0 0];
    layup = [idx ; tply ; theta];
    disp('>')
    disp('> Laminate selected: GLARE 2A-5/4-0.4')
end

% 7) CLS initial geometry
if clsSelect == 1
    % Paper: An ASTM Round Robin (CLS-A) (Johnson 1986)
    lac0 = 0.305-0.254;
    lcb0 = 0.254;
    d = 0.0254;
    b0 = 0.0;
elseif clsSelect == 2
    % CLS Production series 2
    lac0 = 0.06835;
    lcb0 = 0.145;
    d = 0.04366;
    b0 = 0.0;
elseif clsSelect == 3
    % CLS Production series 4
    lac0 = 0.06835+0.110;
    lcb0 = 0.110;
    d = 0.04366;
    b0 = 0.03;
end

% 5) Fatigue load cycle
if lcSelect == 1
    % BOPACS CLS Specimen, K. Hoidus [N]
    Load = LoadCycle('N', 60e3, 0.1, d, sum(tply));
elseif lcSelect == 2
    % Thesis  WSLS Specimen, R. Hanx [N/m]
    Load = LoadCycle('N/m', 230e3, 0.1, d, sum(tply));
elseif lcSelect == 3
    % An ASTM Round Robin CLS specimen, Johnson (1986) [N]
    Load = LoadCycle('N/m', 4.348e3, 0.1, d, sum(tply));
elseif lcSelect == 4
    % Thesis CLS, I. van Teeseling [N]
    Load = LoadCycle('N', 26e3, 0.1, d, sum(tply));
end
disp('>')
disp(['> Load cycle created: ' num2str(Load.P(2)) '[N] (min) - ' num2str(Load.P(1)) '[N] (max)'])

%% Model Code - DO NOT CHANGE!

% Create the composite laminate
materials = [[alu.E1 alu.E2 alu.G alu.v12 alu.ct1 alu.ct2]' ...
    [gf.E1 gf.E2 gf.G gf.v12 gf.ct1 gf.ct2]'];

% Adherent (ad) laminate
adFML = Laminate(materials,  ...
    layup,  ...
    'eMethod', 'membrane', ...
    'sState', 'Plane Stress',  ...
    'dT', 0);

% Adherent meshing (xAC, xCB, xAB and xBC)
adMesh = clsMeshing(lac0,  ...
    lcb0,  ...
    [mesh.qac mesh.qcb], ...
    'qcrack', mesh.qcrack,  ...
    'pasval', NaN,  ...
    'xloc', mesh.xpos);

% Adherent stiffness
adStiff = AdherentStiffness(adFML.Ex,  ...
    adh.E1,  ...
    adFML.h,  ...
    0);

% Overlap (ol) edge loads
olEdgeLoads = clsFBDSolver(adMesh.xAC, ...
    adMesh.xCB,  ...
    Load.Prunning, ...
    adStiff.EIxxac,  ...
    adStiff.EIxxcb,  ...
    adMesh.lAC,  ...
    adMesh.lCB, ...
    adFML.h,  ...
    0,  ...
    bcSelect);

% Overlap adhesive stresses and adherent loads
% Note: The solution of Luo & Tong utilizes the xBC reference frame and the
%       applied load F which is the horizontal load component of P where P
%       is applied along the neutral line of the undeformed CLS

% Get the horizontal load P component; F
alpha = (adFML.h+0)./(2*(lac0+lcb0));
F = Load.Prunning*cos(alpha);

% Overlap region top and bottem adherent load distriubtions
olLoadDistr = clsOverlapLoads(adMesh.xBC, ...
    F, ...
    olEdgeLoads.Mk, ...
    olEdgeLoads.Vk, ...
    adMesh.lCB, ...
    adFML.Ex, ...
    adFML.h, ...
    adh.E1, ...
    adh.G, ...
    0.133e-3);

% Mechanical strain - AC section (out = eng. strain)
em_ac = MechStrainCycle(adFML.plyZ, ...
    adStiff.EAxxac, ...
    adStiff.EIxxac, ...
    Load.Prunning, ...
    olEdgeLoads.Mac, ...
    [1 size(layup, 2)]);

% Mechanical strain - CB section; top adherent (out = eng. strain)
em_cbt = MechStrainCycle(adFML.plyZ, ...
    adStiff.EAxxac, ...
    adStiff.EIxxac, ...
    olLoadDistr.Nt, ...
    olLoadDistr.Mt, ...
    [1 size(layup, 2)]);

% Mechanical strain - CB section; bottom adherent (out = eng. strain)
em_cbb = MechStrainCycle(adFML.plyZ, ...
    adStiff.EAxxac, ...
    adStiff.EIxxac, ...
    olLoadDistr.Nb, ...
    olLoadDistr.Mb, ...
    [1 size(layup, 2)]);

% Total cycle strain - AC section (out = true strain)
%ethreshold = adFML.epth(1,1,2)-adFML.epth(1,1,1);
ethreshold = Inf;
e_ac = MetalStrainCycle([0 0 ; 4.54e-3 325.98e6 ; 9.87e-3 364.91e6], ...
    em_ac.exx, adFML.epres(1,1,1), ...
    ethreshold, ...
    'eMethod', 'TrueStressStrain', ...
    'eTypeInp', 'EngStrain', ...
    'eTypeOut', 'EngStrain');

% Total cycle strain - CB section; top adherent (out = true strain)
e_cbt = MetalStrainCycle([0 0 ; 4.54e-3 325.98e6 ; 9.87e-3 364.91e6], ...
    em_cbt.exx, adFML.epres(1,1,1), ...
    ethreshold, ...
    'eMethod', 'TrueStressStrain', ...
    'eTypeInp', 'EngStrain', ...
    'eTypeOut', 'EngStrain');

% Total cycle strain - CB section; bottom adherent (out = true strain)
e_cbb = MetalStrainCycle([0 0 ; 4.54e-3 325.98e6 ; 9.87e-3 364.91e6], ...
    em_cbb.exx, adFML.epres(1,1,1), ...
    ethreshold, ...
    'eMethod', 'TrueStressStrain', ...
    'eTypeInp', 'EngStrain', ...
    'eTypeOut', 'EngStrain');

% Metal ply stress cycle - AC section (out = true stress)
S_ac = MetalStressCycle([0 0 ; 4.54e-3 325.98e6 ; 9.87e-3 364.91e6], ...
    e_ac.exx, ...
    'eTypeIn', 'TrueStrain', ...
    'sTypeOut', 'EngStress');

% Metal ply stress cycle - CB section; Top adherent (out = true stress)
S_cbt = MetalStressCycle([0 0 ; 4.54e-3 325.98e6 ; 9.87e-3 364.91e6], ...
    e_cbt.exx, ...
    'eTypeIn', 'TrueStrain', ...
    'sTypeOut', 'EngStress');

% Metal ply stress cycle - CB section; Bottom adherent (out = true stress)
S_cbb = MetalStressCycle([0 0 ; 4.54e-3 325.98e6 ; 9.87e-3 364.91e6], ...
    e_cbb.exx, ...
    'eTypeIn', 'TrueStrain', ...
    'sTypeOut', 'EngStress');

% Load DAF % SERR Footprint
load('SERR_DIFF_Fc130.mat');

serr_diff_Fc130(:,1,:) = (serr_diff_Fc130(:,1,:)-55)/1000;
a = serr_diff_Fc130(:,3,:)/100;
b = serr_diff_Fc130(:,4,:)/100;
serr_diff_Fc130(:,2,:) = b;
serr_diff_Fc130(:,3,:) = a;

% Strain Energy Release Rate - DAF effect optionally included
serr = SERRCalculator('Fern1und1991', ...
    'P', Load.Prunning, ...
    'Mk', olEdgeLoads.Mk, ...
    'Mk0', olEdgeLoads.Mk0, ...
    'EAxxac', adStiff.EAxxac, ...
    'EIxxac', adStiff.EIxxac, ...
    'xserr', adMesh.xCB(1,1:mesh.qcrack+1)', ...
    'DAF', serr_diff_Fc130(:,:,2), ...
    'x0daf', 0.005);

% Average disbond growth rate per disbond increment
[dbdN, dG1eq] = adhesiveDGR(serr.GIdaf, ...
    serr.GIIdaf, ...
    adh.c0, ...
    adh.m0, ...
    adh.c100, ...
    1e-10);

% Metal fatigue accumulation per disbond increment
% Note: only fatigue accumulation evaluation in the aluminum cover sheet
%       adjacent to the adhesive bond line of the top (longer) adherent
Saxxab = matrixZipper(S_ac.Sxxa, S_cbt.Sxxa, NaN, 'zip');
Smxxab = matrixZipper(S_ac.Sxxm, S_cbt.Sxxm, NaN, 'zip');
[Minor, dMinor, dN, N, Nf] = metalFatigueInitiation('MH - Sheet', ...
    Saxxab, ...
    Smxxab, ...
    dbdN, ...
    adMesh.xCB(1,2)-adMesh.xCB(1,1));

%% Plotting

gui = 1;

if gui == 1
    % Store in application data to allow acces by the GUI
    setappdata(0, 'l_A', lac0*1000);
    setappdata(0,'N1', olLoadDistr.Nt);
    setappdata(0,'Q1', olLoadDistr.Qt);
    setappdata(0,'M1', olLoadDistr.Mt);
    setappdata(0,'N2', olLoadDistr.Nb);
    setappdata(0,'Q2', olLoadDistr.Qt);
    setappdata(0,'M2', olLoadDistr.Mt);
    setappdata(0,'Shear_a', olLoadDistr.Sxya);
    setappdata(0,'Peel_a', olLoadDistr.Syya);
    setappdata(0,'N', N);
    setappdata(0,'b', adMesh.b*1000);
    setappdata(0,'dbdN', dbdN(1:length(dN)));
    setappdata(0,'dMinor', dMinor);
    setappdata(0,'Minor', Minor);
    setappdata(0,'Sm_nom_ACB', Smxxab);
    setappdata(0,'Sa_nom_ACB', Saxxab);
    setappdata(0,'Sm_nom_BT', S_cbt.Sxxm);
    setappdata(0,'Sa_nom_BT', S_cbt.Sxxa);
    setappdata(0,'Sm_nom_BL', S_cbb.Sxxm);
    setappdata(0,'Sa_nom_BL', S_cbb.Sxxa);
    setappdata(0,'Sm_nom_A', S_ac.Sxxm);
    setappdata(0,'Sa_nom_A', S_ac.Sxxa);
    setappdata(0,'GI', serr.GI);
    setappdata(0,'GII', serr.GII);
    setappdata(0,'dG1_eq', dG1eq);
    setappdata(0,'xAB', adMesh.xAB*1000);
    setappdata(0,'x', adMesh.xBC*1000);
    setappdata(0,'Sy', alu.sy1);
    setappdata(0,'Su', alu.su1);
    
    run('GUIResultsPlotter');
end

% figure(length(findobj('type','figure'))+1)
% hold on
% for n = [1 500 800]
%     plot(adMesh.xAC(n,:), e_ac.exx(n,:,2,1), 'b')
%     plot(adMesh.xAC(n,:), e_ac.exxm(n,:,2,1), 'r')
%     plot(adMesh.xAC(n,:), e_ac.exxpl(n,:,2,1), 'g')
%     plot(adMesh.xAC(n,:), e_ac.exxr(n,:,2,1), 'c')
%     plot(adMesh.xCB(n,:)+adMesh.lAC(n), e_cbt.exx(n,:,2,1), '--b')
%     plot(adMesh.xCB(n,:)+adMesh.lAC(n), e_cbt.exxm(n,:,2,1), '--r')
%     plot(adMesh.xCB(n,:)+adMesh.lAC(n), e_cbt.exxpl(n,:,2,1), '--g')
%     plot(adMesh.xCB(n,:)+adMesh.lAC(n), e_cbt.exxr(n,:,2,1), '--c')
% end
% hold off
%
% figure(length(findobj('type','figure'))+1)
% nn = 1000;
% hold on
% plot(adMesh.xAC(nn,:), S_ac.Sxx(nn,:,2,1), 'b')
% plot(adMesh.xAC(nn,:), S_ac.Sxxm(nn,:,2,1), 'r')
% plot(adMesh.xAC(nn,:), S_ac.Sxxa(nn,:,2,1), 'g')
% plot(adMesh.xCB(nn,:)+adMesh.lAC(nn), S_cbt.Sxx(nn,:,2,1), '--b')
% plot(adMesh.xCB(nn,:)+adMesh.lAC(nn), S_cbt.Sxxm(nn,:,2,1), '--r')
% plot(adMesh.xCB(nn,:)+adMesh.lAC(nn), S_cbt.Sxxa(nn,:,2,1), '--g')
% hold off
