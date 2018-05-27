function [abd, LmM, LmB] = ABD_Matrix_Generator(Al, GF, skn)
%% Description
% Syntax 
%   >ABD_Matrix_Generator(Al, GF, skn)
% Input  
%   >Al     : Structurer with aluminium ply properties
%               .Ex  = Youngs modulus in x direction
%               .Ey  = Youngs modulus in y direction
%               .G   = Shear modulus
%               .v12 = Poisson ratio
%               .ct1 = Thermal coefficient in the x direction
%               .ct2 = Thermal coefficient in the x direction 
%   >GF     : Structurer with glass fibre ply properties
%               .E1 = Youngs modulus in x direction
%               .E2 = Youngs modulus in y direction
%               .G  = Shear modulus
%               .v12 = Poisson ratio
%               .ct1 = Thermal coefficient in the x direction
%               .ct2 = Thermal coefficient in the x direction 
%   >skin   : Structurer the laminate lay-up properties
%               .layer = [1xn] vector with index indicating the ply
%               material; n = # of plies
%                   1 = Aluminium
%                   2 = Glass fibre
%               .t = [1xn] vector with the ply thicknesses; n = # of plies
%               .theta = ; [1xn] vector with the ply angles (ply 1-axis
%               w.r.t. the global, laminate x-axis); n = # of plies
% Output
%   >abd    : Structurer
%               .ABD    = Laminate ABD matrix
%               .comp   = Laminate compliance matrix
%               .stiff  = Laminante ply stiffnesses
%               .zply   = Ply Z-coordinates 
%   >LmM    : Laminate Young's Modulus
%   >LmB    : Laminate Bending Modulus
% Description
%
% -------------------------------------------------------------------------
%
%% Laminate Geometry
% Pre-allocate memory
h = sum(skn.t);
z           = [-h/2 zeros(1,length(skn.theta))];
M           = zeros(3,3,length(skn.theta));
Stiffness   = zeros(3,3,length(skn.theta));
ct          = zeros(3,1,length(skn.theta));

% Create stiffness matrix for each ply
for i = 1:length(skn.theta)
    if skn.layer(i) == 1
        % Aluminium ply
        % Ply stiffness matrix
        [Qmatrix] = QMatrix(Al.E1, Al.E2, Al.G, Al.v12);
        % Ply z-coordinate
        z(i+1)=z(i)+skn.t(i);
        % Transform the thermal exp. coeficients to the global laminate axis
        ct(:,:,i) = skn.t(i)/h*TransformThermalC(Al.ct1, Al.ct2, skn.theta(i));
    elseif skn.layer(i) == 2
        % Glass fibre ply
        % Ply stiffness matrix
        [Qmatrix] = QMatrix(GF.E1, GF.E2, GF.G, GF.v12);
        % Ply z-coordinate
        z(i+1)=z(i)+skn.t(i);
        % Transform the thermal exp. coeficients to the global laminate axis
        ct(:,:,i) = skn.t(i)/h*TransformThermalC(GF.ct1, GF.ct2, skn.theta(i));
    end
    % Transformation matrix
    M(:,:,i) = TransformationMatrix(skn.theta(i));
    
    % Transform local ply stiffness to the global laminate axis
    Stiffness(:,:,i) = M(:,:,i)*Qmatrix*transpose(M(:,:,i));
end

%% Global Stiffness Matrix for Thick Laminates (ABD Matrix)

% Pre-allocate memory
Am = zeros(3,3); 
Bm = zeros(3,3); 
Dm = zeros(3,3);

for i = 1:length(skn.theta)
    % Ply contribution to the A values
    Am = Am + Stiffness(:,:,i)*(z(i+1)-z(i));
    % Ply contribution to the B values
    Bm = Bm + 1/2*Stiffness(:,:,i)*(z(i+1)^2-z(i)^2);
    % Ply contribution to the D values
    Dm = Dm + 1/3*Stiffness(:,:,i)*(z(i+1)^3-z(i)^3);
end

% Construct the ABD matrix
ABD = [Am Bm;Bm Dm];
% Remove singularities
ABD(abs(ABD)<1e-5) =0;
% Compliance matrix
Smat = ABD^(-1);

%% Outputs
% ABD and compliance matrix
abd.ABD     = ABD;   
abd.comp    = Smat;
abd.CT      = ct;
abd.stiff   = Stiffness;
abd.zply    = z;

% Laminate Young's Modulus
LmM.Ex  = 1/(h*Smat(1,1));
LmM.Ey  = 1/(h*Smat(2,2));
LmM.Gxy = 1/(h*Smat(3,3));
LmM.vxy = -Smat(1,2)/Smat(1,1);
LmM.vyx = -Smat(1,2)/Smat(2,2);
% Laminate Bending Modulus
LmB.Ex  = 12/(h^3*Smat(4,4));
LmB.Ey  = 12/(h^3*Smat(5,5));
LmB.Gxy = 12/(h^3*Smat(6,6));
LmB.vxy = -Smat(4,5)/Smat(4,4);
LmB.vyx = -Smat(4,5)/Smat(5,5);
end

function [Q] = QMatrix(E1,E2,G12,v12)
%% Ply stiffness matrix (Thin Classical Laminate Theory)

% Ply properties
v21 = E2/E1*v12;

Q11 = E1/(1-v12*v21);
Q12 = v12*E2/(1-v12*v21);
Q13 = 0;
Q21 = v12*E2/(1-v12*v21);
Q22 = E2/(1-v12*v21);
Q23 = 0;
Q31 = 0;
Q32 = 0;                       
Q33 = G12;                     

% Construct the ply stiffness matrix
Q = [Q11 Q12 Q13;  
     Q21 Q22 Q23;  
     Q31 Q32 Q33];
end

function M = TransformationMatrix(theta)
%% Transformation matrix

% Prepare input
theta = deg2rad(theta);
c = cos(theta);
s = sin(theta);

% Construct transformation matrix
M = [c^2 s^2 2*c*s;
    s^2 c^2 -2*c*s;
    -c*s c*s c^2-s^2];
end

function CT = TransformThermalC(ct1, ct2, theta)
%% Transformation of the Thermal Expansion Coefficients
% Transform from the ply 1,2 coordinate system to the laminate global x,y
% coordinate system

c = cos(theta);
s = sin(theta);

CT = [ct1*c^2+ct2*s^2 ; ct1*s^2+ct2*c^2 ; (ct1-ct2)*c*s];

end
