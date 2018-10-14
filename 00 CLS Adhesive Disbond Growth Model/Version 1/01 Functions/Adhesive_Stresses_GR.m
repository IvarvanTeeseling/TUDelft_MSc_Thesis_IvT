function [Shr_a, Peel_a] = Adhesive_Stresses_GR(x, P, l_B, E, t, E_a, G_a, t_a, v)

x = x./l_B;

P = P/t;

u2c = sqrt(3*(1-v^2)/2)*l_B/t.*sqrt(P/E);

k = 1./(1+2*sqrt(2)*tanh(u2c));

beta = sqrt(8*G_a./E.*t/t_a);

y = sqrt(sqrt(6*E_a/E*t/t_a));

lambda = y*l_B/t;

A = 1/2*(sinh(2*lambda)+sin(2*lambda));

R1 = cosh(lambda).*sin(lambda)+sinh(lambda).*cos(lambda);

R2 = sinh(lambda).*cos(lambda)-cosh(lambda).*sin(lambda);

H = 1./A.*(t./l_B).^2.*P;

Shr_a = (beta.*l_B/t.*(1+3*k).*cosh(l_B.*beta.*x/t)./sinh(beta.*l_B/t)+3*(1-k)).*P.*t./(8*l_B);

kq = k.*l_B/t.*sqrt(3*(1-v^2)*P/E);

Peel_a = H.*((R2.*lambda.^2.*k/2-lambda.*kq.*cosh(lambda).*cos(lambda)).*cosh(lambda.*x).*cos(lambda.*x)+...
    (R1.*lambda.^2.*k/2-lambda.*kq.*sinh(lambda).*sin(lambda)).*sinh(lambda.*x).*sin(lambda.*x));
end