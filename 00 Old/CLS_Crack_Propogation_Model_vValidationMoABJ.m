clear all; close all; clc

%% Input

% Load
S_max = 200e8; % [N/m^2]
R_load = 0.1; % [-]

% Total number of cycles
N = 1e9; % [-]

% Adherent
E   = 70e9; % [Pa]
t   = 0.0016; % [m]
v   = 0.34; % [-]
% Strap (upper, long adherent) length
L_1 = 32*0.0016*(1+1.25); % [m]
% lap (lower, short adherent) length
L_2 = 32*0.0016; % [m]

% Adhesive
t_a = 0.078*0.0016; % [m]
E_a = 0.04*70e9; % [Pa]
G_a = 1.1e9; % [Pa]
v_a = 0.4; % [-]

% Initial crack length
b = 0; % [m]

%% Calculations

% Initital geometry
l_A = L_1-(L_2+b);
l_B = L_2-b;

% Initial angle neutral axis w.r.t. adherent neutral axis
alpha = (t + t_a)/(2*(l_A+l_B));

% Applied force
P_max = S_max*t; % [N/m]
P_min = P*R_load; % [N/m]
P_mean = P*(1+R_load)/2; % [N/m]

P = [P_min P_mean P_max];

% Applied stress
S_min = S_max*R_load; % [N/m^2]
S_mean = S_max*(1+R_load)/2; % [N/m^2]

% Adherent stiffnesses
% Free adherent
D_11 = E*t^3/12;
A_11 = E*t;
% Overlap region (lumped together)
D_00 = 2*E*(t^3/12+t*(t/2+t_a/2)^2) + E_a*t_a^3/12;
A_00 = E*t*2 + E_a*t_a;

% x-vector spanning the free adherent length 
x1 = 0:l_A/100:l_A;
% x-vector spanning the overlap length
x0 = 0:l_B/100:l_B;
% x-vector spanning the entire length
x = [x1 x0+l_A];

% Calculate displacement and overlap edge loads

% >> Roller-Roller Boundary Conditions

P_in =(8/l_B)^2*D_11;
P = P_in/cos(alpha);

% Eigenvalues
lambda_1 = sqrt(P/D_11);
lambda_0 = sqrt(P/D_00);

% Integration coeficients
A_1 = 0;
B_1 = -sqrt(D_11) * (cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) * t + cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) * t_a + (2 * alpha * l_A) + 0.2e1 * alpha * l_B - t - t_a) / (sqrt(D_11) * cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) * sinh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) + sqrt(D_00) * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B)) / 0.2e1;
A_0 = (-0.2e1 * sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) * D_11 ^ (-0.1e1 / 0.2e1) * sinh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) + (-0.2e1 * l_A - 0.2e1 * l_B) * alpha + t + t_a) / cosh(D_00 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_B) / 0.2e1;
B_0 = sqrt(D_00) * B_1 * cosh(D_11 ^ (-0.1e1 / 0.2e1) * sqrt(P) * l_A) * D_11 ^ (-0.1e1 / 0.2e1);

% From force equilibrium (M_A = M_B = 0) due to roller BC
R_A = 0;
M_A = 0;

% Vertical displacement of the free adherent
w1 = A_1 * cosh(lambda_1.*x1) + B_1 * sinh(lambda_1.*x1) + (alpha + R_A/P) .* x1 + M_A/P;
% Vertical desiplacement of the overlap region
w0 = A_0 * cosh(lambda_0.*x0) + B_0 * sinh(lambda_0.*x0) + alpha .* (l_A + x0) - (t + t_a)/2 + (l_A + x0) * R_A/P + M_A/P;

% > Validation < of the displacement using MoABJ
xvalues = [-0.998230721115182 -0.950661883911817 -0.9005706115883939 -0.8479058427052352 -0.7979473301249729 -0.7454153209849752 -0.700797834994957 -0.6535966401572694 -0.6064260821833872 -0.49907451140584924 -0.40230287092944583 -0.29253098791120435 -0.19848602831356832	-0.10707583900328066 -0.04178868223189536 0.002604133423541821];
yvalues = [ -0.16484356051418875 -0.1258275144567702 -0.09384518171489853 -0.06499036215326093 -0.04316414976320254 -0.024465450553378343 -0.011227772317040474 -0.00033636723387414635 0.008211317768104465 0.02060648224976705 0.023636978694614313 0.02118602959010432 0.015625438809247244 0.008504921046248948 0.0029724140572143787 -0.000977188301824139];

figure(1)
hold on
plot((-l_B:l_B/100:0)/l_B,-w0/t,'r')
plot(xvalues,yvalues,':b')
hold off

% Bending moment of the free adherent
% TO DO: Verifiy the sign (normal it is -, but because we define downwards
% displacement as positive, I think it should be +)
M1 = P * (-A_1 * cosh(lambda_1 .* x1) - B_1 * sinh(lambda_1 .* x1));
% Bending moment of the overlap region
M0 = P * (-A_0 * cosh(lambda_0 .* x0) - B_0 * sinh(lambda_0 .* x0));

% Shear force in the free adherent
Q1 = P * (-A_1*lambda_1*sinh(lambda_1.*x1)-B_1*lambda_1*sinh(lambda_1.*x1));
% Shear force in the overlap region
Q2 = P * (-A_0*lambda_0*sinh(lambda_0.*x0)-B_0*lambda_0*sinh(lambda_0.*x0));

% Overlap edge loads
% > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 31
M_k = M1(end);
Q_k = Q1(end);

% Adhesive Stresses
F=P_in;
x0 = -l_B:l_B/1000:0;

% > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 30
beta_k = sqrt(F/D_11);
beta_0 = sqrt(F/D_00);

% > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 33
alpha_k = A_11 * t^2 / (4 * D_11);
alpha_a = (1 + alpha_k) / 4;

% Bending moment factor
k = M_k/(P_in*(t+t_a)/2);
% By GR (note: t_a = 0 in for this expression)
% > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 40
k = 1/(1+2*sqrt(2)*tanh(beta_k*l_B/(2*sqrt(2))*coth(beta_k*l_A)));

% Peel stress integration constants 
% > From Luo and Tong 2004
beta_s = sqrt(2)/2*(24*E_a/(E*t^3*t_a))^(1/4);
beta_t = sqrt(8*G_a/(A_11*t_a));
beta_a = sqrt(alpha_a*beta_t^2);

% Peel stress integration constants 
% > From Luo and Tong 2004
B_s1 = 6/(beta_s^3*E*t^3)*(M_k*beta_s*(sinh(beta_s*l_B)*cos(beta_s*l_B)+cosh(beta_s*l_B)*sin(beta_s*l_B))+Q_k*sinh(beta_s*l_B)*sin(beta_s*l_B))/(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));
B_s4 = 6/(beta_s^3*E*t^3)*(M_k*beta_s*(sinh(beta_s*l_B)*cos(beta_s*l_B)-cosh(beta_s*l_B)*sin(beta_s*l_B))+Q_k*cosh(beta_s*l_B)*cos(beta_s*l_B))/(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));

% Adhesive shear stress (Luo and Tong 2007)
Shear_a1 = beta_t*(F*t+6*M_k)*cosh(beta_t*x0)/(8*t*sinh(beta_t*l_B)) + 3*(F*t-2*M_k)/(8*t*l_B);
% Adheisve peel stress (Luo and Tong 2004)
Peel_a1 = 2*E_a/t_a*(B_s1*sinh(beta_s*x0).*sin(beta_s*x0)+B_s4*cosh(beta_s*x0).*cos(beta_s*x0));

G_I_max = t_a/(2*E_a)*Peel_a1(1)^2;
G_II_max = t_a/(2*G_a)*Shear_a1(1)^2;
G_t_max = G_I_max+G_II_max;
MR = G_II_max/G_t_max

G_eq1_max = sqrt(G_I_max)/2+sqrt(G_I_max/2+G_II_max);
G_eq1_min = sqrt(G_I_min)/2+sqrt(G_I_min/2+G_II_min);
dG_eq1 = (G_eq1_max-G_eq1_min)^2;

c0 = 5.27*10^(-17);
m0 = 3.78;
c100 = 2.23*10^(-24);

k1 = (1+t_a/t)*k;
Shear_a_max = 1/8*((3*k1+1)*beta_t*coth(beta_t*l_B)+3*(1-k1))*F;
Peel_a_max = (1+(beta_k/beta_s)*coth(beta_k*l_A))*beta_s^2*M_k;



% > Validation < of the adhesive stresses using MoABJ
values= [-0.8006306983706201;4.064671588496779
-0.8110654444860231;4.049983038221711
-0.8210032979292641;4.03599394272166
-0.8391398804631789;4.010463843434053
-0.8470907227815917;4.224497004584933
-0.8550410055361846;4.213305728184892
-0.8632408537544983;4.6522135994992055
-0.8711916960729111;4.866246760650085
-0.879142538391324;5.080279921800965
-0.8851058100210886;5.297110902051855
-0.8908211948785921;5.739516047241182
-0.8970340319720778;6.40644617520644
-0.9027499763934014;7.074075757946687
-0.908962253923067;7.515781448361011
-0.9149272042442915;8.408285741264706
-0.9211417200292372;9.750889181882727
-0.9251196592256337;10.871415731437352
-0.928848592522129;11.767067570828544
-0.9333234243906874;12.886894665608168
-0.9368055900425618;14.458569545039651
-0.9412821006025803;16.254069952472037
-0.9450132721543556;18.050619542066954
-0.9489940091698521;20.297268279376226
-0.9529758653129885;22.994365891787353
-0.956709834683864;25.917037669136903
-0.9611908217544424;29.51433357697674
-0.9651765948443188;33.78800225224437
-0.9691634870618353;38.51211980261385
-0.9729036116347307;43.91226039302363
-0.97514690298912;46.8370305346982
-0.9771428671350683;50.21259927886213
-0.9813843608906857;57.413835914904325
-0.9851317597932412;65.74189419347618
-0.9891315219786176;75.64617380751704
-0.9931352011107342;87.12702448441439
-0.9973918030894917;100.40932093433167
-0.9989015062758579;108.06485344673825
-1.0001526909773832;111.66669581061561];

xvalues = values(1:2:end);
yvalues = values(2:2:end);

figure(2)
hold on
plot(x0/l_B,Shear_a1/10^6,'b')
plot(xvalues,yvalues,':r')
hold off

figure(3)
hold on
plot(x0/l_B,Peel_a1/10^6,'b')
hold off

% >> Clamped-Clamped Boundary Conditions