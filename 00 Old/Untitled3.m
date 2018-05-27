clear all;
close all;
clc;

E   = 72e9;
t   = 0.0016; % [m]
v   = 0.33; % [-]

L_1 = 0.2;
L_2 = 0.00635;
b0  = 0;

ta  = 0.15e-3;  % [m]
Ea  = 1.9e9; % [Pa]
Ga  = 0.7e9; % [Pa]
va  = 0.4;

% Initital geometry
l_A = L_1-(L_2-b0);
l_B = L_2-b0;

P = 7500/(0.0254*t);

%%

u2c = sqrt(3*(1-v^2)/2)*l_B/t*sqrt((P)/E)
k = 1/(1+2*sqrt(2)*tanh(u2c))

beta = sqrt(8*Ga*t/(E*ta))

y = sqrt(sqrt(6*Ea*t/(E*ta)))

lambda = y*l_B/t;

A = 0.5*sinh(2*lambda)+sin(2*lambda);

R1 = cosh(lambda)*sin(lambda)+sinh(lambda)*cos(lambda)

R2 = sinh(lambda)*cos(lambda)-cosh(lambda)*sin(lambda)

H = 1/A * (t/l_B)^2*P

m = -1:0.01:0;

shear_a = (beta*l_B/t*(1+3*k)*cosh(l_B*beta*m/t)./sinh(beta*l_B/t)+3*(1-k))*P*t/(8*l_B);

figure(1)
plot(m, shear_a)

kq = k*l_B/t*sqrt(3*(1-v^2)*P/E)

peel_a = H*((R2*lambda^2*k/2-lambda*kq*cosh(lambda)*cos(lambda))*cosh(lambda*m).*cos(lambda*m)+...
    (R1*lambda^2*k/2-lambda*kq*sinh(lambda)*sin(lambda))*sinh(lambda*m).*sin(lambda*m));

%%

P = P*t;

u1 = sqrt(12*(1-v^2)*P/(E*t^3));
u2 = u1/(2*sqrt(2));
k = 1/(1+2*sqrt(2)*tanh(u2*l_B));


