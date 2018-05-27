function [ABD] = ABD_Matrix_Generator()
%% Evalutating Laminate Geometry
h = sum(skn.t);
for i = 1:length(skn.theta)
    z(i) = sum(skn.t(1:i)) - h/2;
end
z = [-h/2 z];

for i = 1:length(skn.theta)
    if (skn.theta(i) == 0) || (skn.theta(i) == 90)
        [Qmatrix] = QMatrix(ud.Ex, ud.Ey, ud.G, ud.vxy);
    else
        [Qmatrix] = QMatrix(pw.Ex, pw.Ey, pw.G, pw.vxy);
    end
    M(:,:,i) = TransformationMatrix(skn.theta(i));
    Stiffness(:,:,i) = M(:,:,i)*Qmatrix*transpose(M(:,:,i));
end

%% Global Stiffness Matrix for Thick Laminates (ABD)
Am = zeros(3,3); Bm = zeros(3,3); Dm = zeros(3,3);

for i = 1:length(skn.theta)
    Am = Am + Stiffness(:,:,i)*(z(i+1)-z(i));
    
    Bm = Bm + 1/2*Stiffness(:,:,i)*(z(i+1)^2-z(i)^2);
    
    Dm = Dm + 1/3*Stiffness(:,:,i)*(z(i+1)^3-z(i)^3);
end

ABDa = [Am Bm;Bm Dm];
ABDa(abs(ABDa)<1e-5) =0;
ABD = ABDa;

vyx = ABD(1,2)/ABD(1,1);
vxy = ABD(1,2)/ABD(2,2);

end

function [Q] = QMatrix(E1,E2,G12,v12)

v21 = E2/E1*v12;

Q = [E1/(1-v12*v21) v12*E2/(1-v12*v21) 0;...
    v12*E2/(1-v12*v21) E2/(1-v12*v21) 0;...
    0 0 G12];

% Compliance = inv(Q);

end

function M = TransformationMatrix(theta)

theta = deg2rad(theta);

c = cos(theta);
s = sin(theta);

M = [c^2 s^2 2*c*s;...
    s^2 c^2 -2*c*s;...
    -c*s c*s c^2-s^2];

end