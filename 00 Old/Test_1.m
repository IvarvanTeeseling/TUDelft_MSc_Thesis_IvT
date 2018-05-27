clear all; close all; clc

v = 0.3;
P = 100000;
t = 0.001;
E = 900e9;
c = 0.5;

%%

%D = E * t^3 / (12 * (1 - v^2));

D = E * t^3 / 12;

%lambda = sqrt(12*(1-v^2)) * sqrt(P / (t*E)) / t;

lambda = sqrt(P/D);

K1 = 1 / (1 + 2 * sqrt(2) * tanh(lambda * c / (2 * sqrt(2))))

%%

u2 = sqrt(3*(1-v^2)/2) * 1/t * sqrt(P/(t*E));

K11 = cosh(u2 * c) / (cosh(u2 * c) + 2 * sqrt(2) * sinh(u2 * c))
