clear all; close all; clc

%% Input
% % Loads (per unit width)
% Load.P   = 200000; 
% 
% % Strap adherent
% Strap.E  = 70e9;     
% % Strap adherent thickness [m]
% Strap.t  = 0.0016;
% % Strap adherent poisson [-]
% Strap.v  = 0.34; 
% % Strap length [m]
% Strap.L  = (32*1.25+32)*0.0016;
% 
% % Lap adherent modulus [Pa]
% Lap.E  = 70e9; 
% % Lap adherent thickness [m]
% Lap.t  = 0.0016;        
% % Lap adherent poisson [-]
% Lap.v  = 0.34;         
% % Lap length [m]
% Lap.L  = 32*0.0016;
% 
% % Adhesive modulus [Pa]
% Adhesive.E  = 0.04*70e9; 
% % Adhesive shear modulus [Pa]
% Adhesive.G  = 0.5*Adhesive.E;
% % Adhesive  thickness [m]
% Adhesive.t  = 0.078*0.0016;        
% % Adhesive poisson [-]
% Adhesive.v  = 0.4;  
% % Initial crack length [m]
% Adhesive.la0 = 0;

% Loads (per unit width)
Load.P   = 500; 

% Strap adherent
Strap.E  = 72e9;     
% Strap adherent thickness [m]
Strap.t  = 0.0032;
% Strap adherent poisson [-]
Strap.v  = 0.33; 
% Strap length [m]
Strap.L  = 0.25;

% Lap adherent modulus [Pa]
Lap.E  = 72e9; 
% Lap adherent thickness [m]
Lap.t  = 0.0032;        
% Lap adherent poisson [-]
Lap.v  = 0.33;         
% Lap length [m]
Lap.L  = 0.20;

% Adhesive modulus [Pa]
Adhesive.E  = 3.1e9; 
% Adhesive shear modulus [Pa]
Adhesive.G  = 1.1e9; 
% Adhesive  thickness [m]
Adhesive.t  = 0.0003;  
% Adhesive poisson [-]
Adhesive.v  = 0.4;  
% Initial crack length [m]
Adhesive.la0 = 0;

% Stress conditions
Conditions.Stress = 0;      % 1 = plane strain, 0 = plane stress
Conditions.BCs    = 1;      % 1 = r-r, 2 = r-c, 3 = c-c.
                            % Note: left is free adherent, right is overlap
if Conditions.Stress == 1
   Strap.E  = Strap.E / (1-Strap.v^2);
   Lap.E    = Lap.E / (1-Lap.v^2);
end  

%% Solutions
%[w1, w0, x] = CLS_Analysis_Lai_et_al_Original(Load, Strap, Lap, Adhesive, Conditions);

[w1, w0, M1, M0, x1, shear1, peel1] = CLS_Analysis_MoABJ(Load, Strap, Lap, Adhesive, Conditions);

%[w11, w00, M11, M00, x11] = CLS_Analysis_Lai_et_al_Adjusted(Load, Strap, Lap, Adhesive, Conditions);

