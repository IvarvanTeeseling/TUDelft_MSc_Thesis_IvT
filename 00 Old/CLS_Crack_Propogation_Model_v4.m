clear all; close all; clc

%% Input

% Total number of cycles
N = 10; % [-]

% Adherent
E   = 70e9; % [Pa]
t   = 0.0032; % [m]
v   = 0.34; % [-]
% Strap (upper, long adherent) length
L_1 = 0.25; % [m]
% lap (lower, short adherent) length
L_2 = 0.2; % [m]

% Adhesive
t_a = 0.0003; % [m]
E_a = 3.1e9; % [Pa]
G_a = 1.1e9; % [Pa]
v_a = 0.4; % [-]

% Load
S_max = 250e3/t; % [N/m^2]
R_load = 0.1; % [-]

% Initial crack length
b = 0; % [m]

plots = 0;

%% Calculations

% Applied force
P_max   = S_max*t;          % [N/m]
P_min   = P_max*R_load;     % [N/m]
P       = [P_min ; P_max];  % [N/m]

% Applied stress
S_min   = S_max*R_load;     % [N/m^2]
S       = [S_min ; S_max];  % [N/m^2]

% Adherent stiffnesses
% Free adherent
D_11 = E*t^3/12;
A_11 = E*t;
% Overlap region (lumped together)
D_00 = 2*E*(t^3/12+t*(t/2+t_a/2)^2) + E_a*t_a^3/12;
A_00 = E*t*2 + E_a*t_a;

N = 10^6;
%b_in = L_2*7/8;
%b_in = 0:b_in/5000:b_in; 

b_hist = zeros(1,N);
 l_B_hist = zeros(1,N);
 G_hist = zeros(3,N);
%G_hist = zeros(3,length(b_in));
%MR_hist = zeros(1,length(b_in));
cnt = 0;
for n = 1:1:N
%for b = b_in 
    
    cnt = cnt + 1;
    
    % Initital geometry
    l_A = L_1-(L_2-b);
    l_B = L_2-b;
    
    % Initial angle neutral axis w.r.t. adherent neutral axis
    alpha = (t + t_a)/(2*(l_A+l_B));
    
    % x-vector spanning the free adherent length
    x1 = (0:l_A/5000:l_A).*ones(2,1);
    % x-vector spanning the overlap length
    x0 = (0:l_B/5000:l_B).*ones(2,1);
    % x-vector spanning the entire length
    x = [x1 x0+l_A].*ones(2,1);
    
    % >> Roller-Roller Boundary Conditions
    
    % Solution eigenvalues
    lambda_1 = sqrt(P/D_11);
    lambda_0 = sqrt(P/D_00);
    
    % Integration contants
    A_1 = [0 ; 0];
    B_1 = -(sqrt(D_11)*(cosh(lambda_0*l_B)*t+cosh(lambda_0*l_B)*t_a+2*alpha*(l_A+l_B)-t-t_a))...
        ./(2*(sqrt(D_11)*cosh(lambda_0*l_B).*sinh(lambda_1*l_A)+sqrt(D_00).*cosh(lambda_1*l_A).*sinh(lambda_0*l_B)));
    A_0 = (-2*sqrt(D_00/D_11)*B_1.*cosh(lambda_1*l_A).*sinh(lambda_0*l_B)-2*alpha*(l_A+l_B)+t+t_a)./(2*cosh(lambda_0*l_B));
    B_0 = sqrt(D_00/D_11)*B_1.*cosh(lambda_1*l_A);
    
    % Support reaction forces (see FBD)
    R_A = 0;
    M_A = 0;
    
    % Vertical displacement
    w1 = A_1.*cosh(lambda_1.*x1)+B_1.*sinh(lambda_1.*x1)+(alpha+R_A./P).*x1+M_A./P;
    w0 = A_0.*cosh(lambda_0.*x0)+B_0.*sinh(lambda_0.*x0)+alpha.*(l_A+x0)-(t+t_a)/2+(l_A + x0).*R_A./P+M_A./P;
    
    % Bending moments
    M1 = P.*(-A_1.*cosh(lambda_1.*x1)-B_1.*sinh(lambda_1.*x1));
    M0 = P.*(-A_0.*cosh(lambda_0.*x0)-B_0.*sinh(lambda_0.*x0));
    
    % Shear forces
    Q1 = P.*(-A_1.*lambda_1.*sinh(lambda_1.*x1)-B_1.*lambda_1.*sinh(lambda_1.*x1));
    Q2 = P.*(-A_0.*lambda_0.*sinh(lambda_0.*x0)-B_0.*lambda_0.*sinh(lambda_0.*x0));
    
    % if plots == 1
    %     figure(1)
    %     hold on
    %     plot(x1, -w1, 'r')
    %     plot(x0+l_A, -w0, 'b')
    %     hold off
    % end
    
    % Overlap edge loads
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 31
    M_k = M1(:,end);
    Q_k = Q1(:,end);
    
    % Adhesive Stresses
    x0 = -l_B:l_B/5000:0;
    F = P*cos(alpha);
    
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 30
    beta_k = sqrt(F/D_11);
    beta_0 = sqrt(F/D_00);
    
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 33
    alpha_k = A_11*t^2/(4*D_11);
    alpha_a = (1+alpha_k)/4;
    
    % Bending moment factor
    k = M_k./(P.*(t+t_a)/2);
    % By GR (note: t_a = 0 in for this expression)
    % > From MoABJ, chapter Analysis of Cracked Lap Shear Joints, page 40
    k = 1/(1+2*sqrt(2)*tanh(beta_k*l_B/(2*sqrt(2)).*coth(beta_k*l_A)));
    
    % Peel stress integration constants
    % > From Luo and Tong 2004
    beta_s = sqrt(2)/2*(24*E_a/(E*t^3*t_a))^(1/4);
    beta_t = sqrt(8*G_a/(A_11*t_a));
    beta_a = sqrt(alpha_a*beta_t^2);
    
    % Peel stress integration constants
    % > From Luo and Tong 2004
    B_s1 = 6/(beta_s^3*E*t^3)*(M_k*beta_s*(sinh(beta_s*l_B)*cos(beta_s*l_B)+cosh(beta_s*l_B)*sin(beta_s*l_B))+Q_k*sinh(beta_s*l_B)*sin(beta_s*l_B))...
        /(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));
    B_s4 = 6/(beta_s^3*E*t^3)*(M_k*beta_s*(sinh(beta_s*l_B)*cos(beta_s*l_B)-cosh(beta_s*l_B)*sin(beta_s*l_B))+Q_k*cosh(beta_s*l_B)*cos(beta_s*l_B))...
        /(sinh(2*beta_s*l_B)+sin(2*beta_s*l_B));
    
    % Adhesive shear stress (Luo and Tong 2007)
    Shear_a = beta_t*(F*t+6*M_k)*cosh(beta_t*x0)/(8*t*sinh(beta_t*l_B)) + 3*(F*t-2*M_k)/(8*t*l_B);
    % Adheisve peel stress (Luo and Tong 2004)
    Peel_a = 2*E_a/t_a*(B_s1*sinh(beta_s*x0).*sin(beta_s*x0)+B_s4*cosh(beta_s*x0).*cos(beta_s*x0));
    
    % SERR
    G_I = t_a/(2*E_a)*Peel_a(:,1).^2;
    G_II = t_a/(2*G_a)*Shear_a(:,1).^2;
    G_T = G_I+G_II;
    
    % Mode ratio
    MR = G_II./G_T;
    
    % Disbond growth rate based an D. Burger (2005)
    % Equivalent mode I
    G1_eq = sqrt(G_I)/2+sqrt(G_I/2+G_II);
    dG1_eq = (G1_eq(2)-G1_eq(1)).^2;
    
    % Paris Law fit coefficients (D. Burger data on FM94)
    % TO DO: (1) include threshold and (2) G_c
    c0 = 5.27*10^(-17);
    m0 = 3.78;
    c100 = 10^(-17.6);
    
    % Note: MR different at P_min and _P_max due to geometric non linearity
    dbdN = c100^MR(2)*c0^(1-MR(2))*dG1_eq.^m0;
    
     b = b + dbdN;
    
     b_hist(n) = b;
%     l_B_hist(n) = l_A;
%     G_hist(:,cnt) = [G_I(2) ; G_II(2); G_T(2)];
%     MR_hist(cnt) = MR(2);
%     
%     dbdN_hist(:,cnt) = dbdN;
%     Shear_hist(:,cnt) = Shear_a(:,1);
    
    if b > L_2/8
        b_hist(:,n+1:end) = [];
        l_B_hist(:,n+1:end) = [];
        G_hist(:,n+1:end) = [];
        N_total = 1:n;
        break    
    end
    
end


% vali = [200.09141975988544; 7.453652089114797
% 190.58736456953733; 10.577236667839756
% 180.3155768202687; 11.27485298423862
% 170.29894771907297; 13.484186235187622
% 160.02683759627553; 14.483114343204235
% 150.26859087844102; 15.19104661252581
% 200.13510137304257; 46.62590432490954
% 190.92706567554717; 53.06993625063467
% 179.28261144049924; 56.75316127541714
% 169.90661913448082; 60.18063663399052
% 158.43431236383077; 62.963364934893974
% 150.34346213051293; 65.21138300930315
% 200.29822237863308; 54.16213776632756
% 191.04577972754063; 62.111868987397884
% 179.39842413073322; 68.50690013674011
% 169.84963961325923; 73.437495802428
% 158.8475952277972; 76.68164808093724
% 150.15761379113957; 78.91763087693658];
% vali_x = (vali(1:2:end)-200)'/1000*-1;
% vali_y = vali(2:2:end)';
% 
% figure(6)
% hold on
% plot(b_in, G_hist)
% plot(vali_x(1:6),vali_y(1:6),'-r')
% plot(vali_x(7:12),vali_y(7:12),'-b')
% plot(vali_x(13:18),vali_y(13:18),'-c')
% hold off

% figure(7)
% plot(b_in, MR_hist)
% 
% figure(9)
% plot(b_in, dbdN_hist)

% figure(4)
% plot(N_total, b_hist)
% 
% figure(5)
% plot(N_total, l_B_hist)
% 
% figure(6)
% plot(N_total, G_hist)

if plots == 1
    figure(2)
    hold on
    plot(x0/l_B,Shear_a/10^6,'b')
    hold off
end

if plots == 1
    figure(3)
    hold on
    plot(x0/l_B,Peel_a/10^6,'b')
    hold off
end

% k1 = (1+t_a/t)*k;
% Shear_a_max = 1/8*((3*k1+1)*beta_t*coth(beta_t*l_B)+3*(1-k1))*F;
% Peel_a_max = (1+(beta_k/beta_s)*coth(beta_k*l_A))*beta_s^2*M_k;

% >> Clamped-Clamped Boundary Conditions