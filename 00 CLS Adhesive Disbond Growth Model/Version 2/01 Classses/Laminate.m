classdef Laminate
    % Laminate Summary of this class goes here
    % Detailed explanation goes here
    
    % >> Source Bibliography
    % >
    % @article{Spronk2015,
    %   title={Predicting fatigue crack initiation in fibre metal laminates
    %   based on metal fatigue test data},
    %   author={Spronk, SWF and {\c{S}}en, I and Alderliesten, RC},
    %   journal={International Journal of Fatigue},
    %   volume={70},
    %   pages={428--439},
    %   year={2015},
    %   publisher={Elsevier}
    % }
    % >
    % @article{Tsai1980,
    %   title={Introduction to composite materials},
    %   author={Tsai, S. and Hahn, H.},
    %   year={1980},
    %   publisher={Technomic Pubishing Company Inc.}
    % }
    
    properties (GetAccess = public, SetAccess = private)
        plyMaterials    % Ply material index
        plyLayup        % Lay-up
        plyZ            % Ply z-coordinates (bottom to top)
        h               % Laminate thickness
        ABD             % ABD matrix
        abd             % abd matrix
        Qlocal          % Ply Q matrix (1,2 reference frame)
        Qglobal         % Ply Q matrix (x,y reference frame)
        Ex              % Laminate smeared Ex
        Ey              % Laminate smeared Ey
        Gxy             % Laminate smeared Gxy
        vxy             % Laminate smeared vxy
        eMethod         % 'membrane' (default) or 'flexural'
        sState          % 'Plane Stress' (optional) or 'Plane Strain'
        dT              % dt = Tcuring - Tapplication
        Nt              % Thermal line load
        Mt              % Thermal bending moment
        epres           % Ply residual strain (x,y reference frame)
        epth            % Ply thermal strain (x,y reference frame)
        elth0           % Laminate curing strain (x,y reference frame)
        klth            % Laminate curing curvature (x,y reference frame)
    end
    
    methods
        function obj = Laminate(material, layup, varargin)
            % > Laminate() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sState
            expectedeMethod = {'membrane' 'flexural'};
            expectedsState = {'Plane Stress' 'Plane Strain'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validmaterial   = @(x) size(x, 1) == 4 || size(x, 1) == 6;
            validlayup      = @(x) size(x, 1) == 3 && ...
                max(x(1,:)) == size(material, 2);
            valideMethod    = @(x) any(strcmp(x, expectedeMethod));
            validsState     = @(x) any(strcmp(x, expectedsState));
            validdT         = @(x) isnumeric(x);
            
            % Add required arguments and their validation functions
            addRequired(p, 'material', validmaterial);
            addRequired(p, 'layup', validlayup);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            %       1. methodEng = 'membrane'
            %       2. sState = 'Plane Stress'
            %       3. dT = [];
            addOptional(p, 'eMethod', 'membrane', valideMethod);
            addOptional(p, 'sState', 'Plane Stress', validsState);
            addOptional(p, 'dT', [], validdT);
            
            % Validate the input arguments
            parse(p, material, layup, varargin{:});
            
            % Perform the Classical laminate Theory
            [CLT, eProp] = Laminate.CLT( ...
                p.Results.material, ...
                p.Results.layup, ...
                p.Results.eMethod, ...
                p.Results.sState);
            
            % Get the ply residual (thermal) strains
            if size(p.Results.material, 1) == 6 && ~isempty(p.Results.dT)
                % Residual thermal line loads
                [Nt, Mt] = Laminate.thermalLineLoads(...
                    CLT.Qglobal, ...
                    p.Results.material(5:6,:), ...
                    p.Results.layup, ...
                    p.Results.dT, ...
                    CLT.zply);
                
                % Residual thermal strains (ply and laminate)
                [epres, epth, elth0, klth] = Laminate.residualStrain( ...
                    p.Results.material(5:6,:), ... 
                    p.Results.layup, ...
                    p.Results.dT, ... 
                    CLT.zply, ... 
                    Nt, ... 
                    Mt, ...
                    CLT.abd);
            else
                if size(material, 1) ~= 6 && ~isempty(dT)
                    warning('dT given, but no material ct properties!')
                end
            end
            
            % Set object properties
            obj.plyMaterials    = p.Results.material;
            obj.plyLayup        = p.Results.layup;
            obj.plyZ            = CLT.zply;
            obj.h               = sum(p.Results.layup(2, :));
            obj.ABD             = CLT.ABD;
            obj.abd             = CLT.abd;
            obj.Qlocal          = CLT.Qlocal;
            obj.Qglobal         = CLT.Qglobal;
            obj.eMethod         = p.Results.eMethod;
            obj.sState          = p.Results.sState;
            obj.Ex              = eProp.Ex;
            obj.Ey              = eProp.Ey;
            obj.Gxy             = eProp.Gxy;
            obj.vxy             = eProp.vxy;
            obj.dT              = p.Results.dT;
            obj.Nt              = Nt;
            obj.Mt              = Mt;
            obj.epres           = epres;
            obj.epth            = epth;
            obj.elth0           = elth0;
            obj.klth            = klth;
        end
    end
    
    methods (Static)
        function [CLT, eProp] = CLT(material, layup, eMethod, sState)
            % > CLT() performs the Classical Laminate Theory (CLT)
            %
            % > Several functions are called that are defined in this
            % static method domain of this class
            %------------------------------------------------------------
            
            % Get the ply local and global stiffness matrices
            [Qlocal, Qglobal] = Laminate.QlocalQglobal(material, layup);
            
            % Get the ply bottom and top z-coordinates
            z = Laminate.plyBottomTopZ(layup(2, :));
            
            % Get the laminate ABD and compliance matrix
            [ABD, abd] = Laminate.ABDMatrix(Qglobal, z);
            
            % Get the laminate engineering (smeared) properties
            eProp = Laminate.laminateEngProps(ABD, ...
                abd, sum(layup(2, :)), eMethod, sState);
            
            % Set outputs
            CLT.ABD = ABD;
            CLT.abd = abd;
            CLT.Qglobal = Qglobal;
            CLT.Qlocal = Qlocal;
            CLT.zply = z;
        end
        
        function [ABD, abd] = ABDMatrix(Qglobal, zply)
            % > ABDMatrix() creates and returns the laminate ABD matrix
            %------------------------------------------------------------
            
            % Pre-allocate memory
            A = zeros(3, 3);
            B = zeros(3, 3);
            D = zeros(3, 3);
            
            % Compute the Am, Bm and Dm components of the ABD matrix
            for i = 1:length(zply)-1
                % Ply contribution to the A values
                A = A + Qglobal(:, :, i)*(zply(i+1)-zply(i));
                
                % Ply contribution to the B values
                B = B + 1/2*Qglobal(:, :, i)*(zply(i+1)^2-zply(i)^2);
                
                % Ply contribution to the D values
                D = D + 1/3*Qglobal(:, :, i)*(zply(i+1)^3-zply(i)^3);
            end
            % Construct the ABD matrix
            ABD = [A B ; B D];
            
            % Remove floating point singularities
            ABD(abs(ABD)<1e-5) = 0;
            
            % Create the compliance matrix
            abd = ABD^(-1);
        end
        
        function [z] = plyBottomTopZ(tply)
            % > plyBottomTopZ() returns the ply bottom and top
            % Z-coordinates w.r.t. the laminate neutral axis
            %------------------------------------------------------------
            
            % Store the ply bottom and top Z-coordinate (bottom to top)
            z = [-sum(tply)/2 cumsum(tply)-sum(tply)/2];
        end
        
        function engProp = laminateEngProps(ABD, abd, h, em, ss)
            % > laminateEngProps() creates and returns the laminate
            % Engineering properties, also commonly known as the 'smeared
            % properties', based on the requested method and stress state
            %------------------------------------------------------------
            
            % Calculate the Engineering properties
            switch em
                case 'membrane'
                    % laminate Young's Modulus
                    if strcmp(ss, 'Plane Stress')
                        % Note: plane stress in the width direction
                        engProp.Ex = 1/(h*abd(1, 1));
                    elseif strcmp(ss, 'Plane Strain')
                        % Note: plane strain in the width direction
                        engProp.Ex = 1/h*ABD(1, 1);
                    end
                    engProp.Ey = 1/(h*abd(2, 2));
                    engProp.Gxy = 1/(h*abd(3, 3));
                    engProp.vxy = -abd(1, 2)/abd(1, 1);
                    engProp.vyx = -abd(1, 2)/abd(2, 2);
                    
                case 'flexural'
                    % laminate Bending Modulus
                    if strcmp(ss, 'Plane Stress')
                        % Note: plane stress in the width direction
                        engProp.Ex = 12/(h^3*abd(4, 4));
                    elseif strcmp(ss, 'Plane Strain')
                        % Note: plane strain in the width direction
                        engProp.Ex = 1/h*ABD(4, 4);
                    end
                    engProp.Ey = 12/(h^3*abd(5, 5));
                    engProp.Gxy = 12/(h^3*abd(6, 6));
                    engProp.vxy = -abd(4, 5)/abd(4, 4);
                    engProp.vyx = -abd(4, 5)/abd(5, 5);
            end
        end
        
        function [Qlocal, Qglobal] = QlocalQglobal(material, layup)
            % > QlocalQglobal() creates and returns the local (ply
            % 1,2 reference frame) and global (laminate x,y reference
            % frame) stiffness matrix of each ply by applying the Thin
            % Classical laminate Theory
            %------------------------------------------------------------
            
            % Extract laminate properties
            idx = layup(1, :);
            theta = layup(3, :);
            
            % Pre-allocate memory
            M = zeros(3, 3, length(theta));
            Qlocal = zeros(3, 3, length(theta));
            Qglobal = zeros(3, 3, length(theta));
            
            % Construct the local and global stiffness matrix for each ply
            for i = 1:length(theta)
                % Extract ply material properties
                E1 = material(1, idx(i));
                E2 = material(2, idx(i));
                G12 = material(3, idx(i));
                v12 = material(4, idx(i));
                
                % Q-matrix elements
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
                
                % Assemble the Q-matrix
                Qlocal(:,:,i) = [Q11 Q12 Q13 ; Q21 Q22 Q23 ; Q31 Q32 Q33];
                
                % Transform the Q-matrix from local (1,2) to global (x,y)
                M(:,:,i) = Laminate.transMatrix(theta(i));
                Qglobal(:,:,i) = M(:,:,i)*Qlocal(:,:,i)* ...
                    transpose(M(:,:,i));
            end
        end
        
        function M = transMatrix(theta)
            % > gettransMatrix() creates and a rotation matrix
            % based on the .... transformation sequence
            %------------------------------------------------------------
            
            % Prepare input
            theta = deg2rad(theta);
            c = cos(theta);
            s = sin(theta);
            
            % Construct transformation matrix
            M = [c^2 s^2 2*c*s;
                s^2 c^2 -2*c*s;
                -c*s c*s c^2-s^2];
        end
        
        function [Nt, Mt] = thermalLineLoads(Qglobal, ct, layup, dT, zply)
            % > thermalLineLoads() find the residual thermal stress
            % resultant (Nt) and the thermal moment (Mt) caused by a dT
            % where dT = Tcuring - Tapplication
            %------------------------------------------------------------
            
            % Extract laminate properties
            idx = layup(1,:);
            theta = layup(3,:);
            
            % Initialize
            Nt = zeros(3,1);
            Mt = zeros(3,1);
            
            % Find the thermal residual line loads
            %   > Source: [Tsai1980; EQ8.48]
            for i = 1:size(layup,2)
                ct1 = ct(1, idx(i));
                ct2 = ct(2, idx(i));
                Q = Qglobal(:,:, idx(i));
                
                % Thermal coefficients form material (1,2) to the lamiante
                % (x,y) reference frame
                ctxy = Laminate.Ct12ToCtxy(ct1, ct2, theta(i));
                
                % Ply thermal strain the laminate (x,y) reference frame
                etxy = ctxy*dT;
                
                % Thermal residual line laods
                Nt = Nt+Q*etxy*(zply(i+1)-zply(i));
                Mt = Mt+1/2*Q*etxy*(zply(i+1)^2-zply(i)^2);
            end
            
            % Remove numerical floating point errors (<1e-10)
            Nt(abs(Nt)<1e-10) = 0;
            Mt(abs(Mt)<1e-10) = 0;
        end
        
        function [epres, epth, elth0, klth] = residualStrain( ...
                ct, layup, dT, zply, Nt, Mt, abdmat)
            % > residualStrain() find the residual thermal strains
            % caused by a dT where dT = Tcuring - Tapplication
            %
            % Output:
            %   - epres   = Ply residual strain ((x,y) RF)
            %               > [epresx epresy epresxy]'
            %   - epth    = Ply thermal strain ((x,y) RF)
            %               > [epthx epthy epthxy]'
            %   - elth0   = Laminate thermal strain at z = 0 ((x,y) RF)
            %               > [elthx elthy elthxy]'
            %   - klth    = Laminate thermal curvature ((x,y) RF)
            %               > [klthx klthy klthxy]'
            %------------------------------------------------------------
            
            % Laminate curing strain and curvature
            elth0 = abdmat(1:3,1:3)*Nt+abdmat(1:3,4:6)*Mt;
            klth = abdmat(4:6,1:3)*Nt+abdmat(4:6,4:6)*Mt;
            
            % Extract laminate properties
            idx = layup(1,:);
            theta = layup(3,:);
            
            % Ply mid z-coordinates
            z = (zply(1:end-1)+zply(2:end))/2;
            
            % Pre-allocate memory
            epth = zeros(3,1,size(layup,2));
            epres = zeros(3,1,size(layup,2));
            
            for i = 1:size(layup,2)
                % Material thermal coefficients
                ct1 = ct(1, idx(i));
                ct2 = ct(2, idx(i));
                
                % Thermal coefficients from material (1,2) to the laminate
                % (x,y) reference frame
                ctxy = Laminate.Ct12ToCtxy(ct1, ct2, theta(i));
                
                % Ply thermal strain ((x,y) reference frame)
                epth(:,:,i) = ctxy*dT;
                
                % Ply residual strain ((x,y) reference frame)
                epres(:,:,i) = elth0+klth*z(i)-epth(:,:,i);
            end
        end
        
        function ctxy = Ct12ToCtxy(ct1, ct2, theta)
            % > Ct12ToCtxy() transforms the material thermal
            % coefficients from the local (1,2) to the laminate global (x,y)
            % reference frame
            %------------------------------------------------------------
            
            % Map thermal coefficient from the material principle (1,2)
            % reference frame to the laminate (x,y) reference frame
            %   > Source: [Spronk2015; EQ1]
            c = cos(theta);
            s = sin(theta);
            ctxy = [ct1*c^2+ct2*s^2 ; ct1*s^2+ct2*c^2 ; (ct1-ct2)*c*s];
        end
        
    end
end

