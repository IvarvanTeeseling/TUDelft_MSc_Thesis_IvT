clear all; close all; clc

E1 = 56106e6;
E2 = 56106e6;
t  = 3.6e-3;
Ea = 3900e6;
Gad = 1300e6;
tad = 0.2e-3;
L  = 0.95;
P  = 420e3;
M  = 10000;
Q  = 0;

k = sqrt(Gad/(E1.*t*tad));
c = (6/(E1*t^3*Ea*tad))^(1/4);

a1 = 0;
a2 = P*k/(2*sinh(k*L/2));
a3 = 0;
a4 = P/2+a1/k*cosh(k*L/2);
a5 = P/2-a1/k*cosh(k*L/2);

b1 = P*t/L*sin(c*L/2)*sinh(c*L/2)/(cos(c*L/2)*sin(c*L/2)+cosh(c*L/2)*sinh(c*L/2));
b2 = P*t/L*cos(c*L/2)*cosh(c*L/2)/(cos(c*L/2)*sin(c*L/2)+cosh(c*L/2)*sinh(c*L/2));
b3 = 0;
b4 = 0;
b5 = Q/2+(b3+b4)/(2*c)*sin(c*L/2)*sinh(c*L/2)+(b4-b3)/(2*c)*cos(c*L/2)*cosh(c*L/2);
b6 = Q/2-(b3+b4)/(2*c)*sin(c*L/2)*sinh(c*L/2)-(b4-b3)/(2*c)*cos(c*L/2)*cosh(c*L/2);

x = -L/2:0.0001:L/2;

M1 = -b2./(2.*c^2).*sin(c.*x).*sinh(c.*x)+b1./(2.*c^2).*cos(c.*x).*cosh(c.*x)-b3./(2.*c^2).*sin(c.*x).*cosh(c.*x)+b4./(2.*c^2).*cos(c.*x).*sinh(c.*x)+b5.*x-t.*a1./(2.*k).*cosh(k.*x)-t.*a2./(2.*k).*sinh(k.*x)-t.*a3./2.*x;
P1 = -a1/k*cosh(k*x)-a2/k*sinh(k*x)+a3*x+a4;

peel = P.*t/L.*(sin(c.*L/2).*sinh(c.*L/2).*sin(c.*x)/(cos(c.*L/2).*sin(c.*L/2)+cosh(c.*L/2).*sinh(c.*L/2))+cos(c.*L/2).*cosh(c.*L/2).*cosh(c.*x).*cosh(c.*x)/(cos(c.*L/2).*sin(c.*L/2)+cosh(c.*L/2).*sinh(c.*L/2)));

figure(1)
plot(x,M1)

figure(2)
plot(x,peel)

figure(3)
plot(x,P1)
