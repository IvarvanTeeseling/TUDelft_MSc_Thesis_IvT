classdef SERRCalculator
    %SERRCalculator Summary of this class goes here
    %   Detailed explanation goes here
    
    % >> Source Bibliography
    % >
    % @article{Fernlund1991,
    %   title={Failure load prediction of structural adhesive joints: Part 1: Analytical method},
    %   author={Fernlund, G and Spelt, JK},
    %   journal={International Journal of Adhesion and Adhesives},
    %   volume={11},
    %   number={4},
    %   pages={213--220},
    %   year={1991},
    %   publisher={Elsevier}
    % }
    % >
    % @article{Edde1992,
    %   title={On the fracture parameters in a clamped cracked lap shear adhesive joint},
    %   author={Edde, F and Verreman, Y},
    %   journal={International journal of adhesion and adhesives},
    %   volume={12},
    %   number={1},
    %   pages={43--48},
    %   year={1992},
    %   publisher={Elsevier}
    % }
    % >
    % @article{Johnson1987,
    %   title={Stress analysis of the cracked-lap-shear specimen: an ASTM round-robin},
    %   author={Johnson, W Steven},
    %   journal={Journal of testing and evaluation},
    %   volume={15},
    %   number={6},
    %   pages={303--324},
    %   year={1987},
    %   publisher={ASTM International}
    % }
    
    properties (GetAccess = public)
        Method      % Calculation method deployed
        GI          % Mode I SERR component excl. DAF effect
        GII         % Mode II SERR component excl. DAF effect
        GT          % Total SERR component excl. DAF effect
        MR          % Mode Ratio (GII/GT) excl. DAF effect
        flagDAF     % DAF flagger (default = 'noDAF')
        GIdaf       % Mode I SERR component incl. DAF effect (optional)
        GIIdaf      % Mode II SERR component incl. DAF effect (optional)
        GTdaf       % Totl SERR component incl. DAF effect (optional)
        MRdaf       % Mode Ratio (GII/GT) incl. DAF effect (optional)
    end
    
    methods
        function obj = SERRCalculator(method, varargin)
            % > SERRCalculator() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            expectedMethod = {'Verreman1992' 'Brussat1977' 'Fern1und1991'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum      	= @(x) isnumeric(x);
            validMethod     = @(x) any(strcmp(x, expectedMethod));
            
            % Add required arguments and their validation functions
            addRequired(p, 'method', validMethod);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            addOptional(p, 'ta', [], validnum);
            addOptional(p, 'Ga', [], validnum);
            addOptional(p, 'Ea', [], validnum);
            addOptional(p, 'Sxya', [], validnum);
            addOptional(p, 'Syya', [], validnum);
            addOptional(p, 'P', [], validnum);
            addOptional(p, 't', [], validnum);
            addOptional(p, 'Mk', [], validnum);
            addOptional(p, 'Mk0', [], validnum);
            addOptional(p, 'EAxxac', [], validnum);
            addOptional(p, 'EIxxac', [], validnum);
            addOptional(p, 'EIxxcb', [], validnum);
            addOptional(p, 'xserr', [], validnum);
            addOptional(p, 'DAF', [], validnum);
            addOptional(p, 'x0daf', [], validnum);
            
            % Validate the input arguments
            parse(p, method, varargin{:});
            
            % Input check spanning more than one paramater
            
            switch method
                case 'Verreman1992'
                    % Input check
                    if any([isempty(p.Results.ta); ...
                            isempty(p.Results.Ga); ...
                            isempty(p.Results.Ea); ....
                            isempty(p.Results.Sxya); ...
                            isempty(p.Results.Syya)])
                        error('Verreman1992 requires the inputs: ta, Ga, Ea, Sxya, Syya')
                    end
                    
                    % SERRCalculator Components
                    [GI, GII, GT, MR] = SERRCalculator.Verreman1992(p.Results.ta, ...
                        p.Results.Ga, ...
                        p.Results.Ea, ....
                        p.Results.Sxya, ...
                        p.Results.Syya);
                    
                case 'Brussat1977'
                    % Input check
                    if any([isempty(p.Results.P); ...
                            isempty(p.Results.t); ...
                            isempty(p.Results.Mk); ...
                            isempty(p.Results.EAxxac); ...
                            isempty(p.Results.EIxxac); ...
                            isempty(p.Results.EIxxcb)])
                        error('Brussat1977 requires the inputs: P, t, Mk, EAxxac, EIxxac, EIxxcb')
                    end
                    
                    % SERRCalculator Components
                    [GI, GII, GT, MR] = SERRCalculator.Brussat1977(p.Results.P, ...
                        p.Results.t, ...
                        p.Results.Mk, ....
                        p.Results.EAxxac, ...
                        p.Results.EIxxac, ...
                        p.Results.EIxxcb);
                    
                case 'Fern1und1991'
                    % Input check
                    if any([isempty(p.Results.P); ...
                            isempty(p.Results.Mk); ...
                            isempty(p.Results.Mk0); ....
                            isempty(p.Results.EAxxac); ...
                            isempty(p.Results.EIxxac)])
                        error('Fern1und1991 requires the inputs: P, Mk, Mk0, EAxxac, EIxxac')
                    end
                    
                    % SERRCalculator Components
                    [GI, GII, GT, MR] = SERRCalculator.Fern1und1991(p.Results.P, ...
                        p.Results.Mk, ...
                        p.Results.Mk0, ....
                        p.Results.EAxxac, ...
                        p.Results.EIxxac);
            end
            
            % Set output
            obj.Method = method;
            obj.GI = GI;
            obj.GII = GII;
            obj.GT = GT;
            obj.MR = MR;
            
            % Incorporate the DAF effect on the SERR (if requested)
            if ~any([isempty(p.Results.xserr), ...
                    isempty(p.Results.DAF), ...
                    isempty(p.Results.x0daf)])
                
                [GIdaf, GIIdaf, GTdaf, MRdaf] = SERRCalculator.SERRDAF( ...
                    p.Results.xserr, ...
                    GI, ...
                    GII, ...
                    p.Results.DAF, ...
                    p.Results.x0daf);
                
                % Set output
                obj.flagDAF = 'onDAF';
                obj.GIdaf = GIdaf;
                obj.GIIdaf = GIIdaf;
                obj.GTdaf = GTdaf;
                obj.MRdaf = MRdaf;
            else
                % Set output
                obj.flagDAF = 'offDAF';
                obj.GIdaf = GI;
                obj.GIIdaf = GII;
                obj.GTdaf = GT;
                obj.MRdaf = MR;
            end
        end
    end
    
    methods (Static)
        function [GI, GII, GT, MR] = Verreman1992(ta, Ga, Ea, Sxya, Syya)
            % > laminateStrainCycle() finds SERRCalculator Mode I and Mode II using
            % the adhesive stress formulation
            %
            % > Source: [Verreman1992]
            %------------------------------------------------------------
            
            % > Source: [Verreman1992]
            GI = ta/(2*Ea)*Sxya.^2;
            
            % > Source: [Verreman1992]
            GII = ta/(2*Ga)*Syya.^2;
            
            GT = GI+GII;
            MR = GII./GT;
        end
        
        function [GI, GII, GT, MR] = Brussat1977(P, t, Mk, EAxxac, ...
                EIxxac, EIxxcb)
            % > Brussat1977() finds SERRCalculator for an infinitely long
            % joint (independend of crack length
            %
            % > Source: [Brussat1977]
            %------------------------------------------------------------
            
            % > Source: [Brussat1977, EQII-1]
            GT = P.^2/(4*EAxxac*t);
            
            % > Source: [Brussat1977, EQII-5]
            GI = 2*Mk.^2/(7*EIxxac)*(1-EIxxac/EIxxcb);
            GII = GT-GI;
            MR = GII./GT;
        end
        
        function [GI, GII, GT, MR] = Fern1und1991(P, Mk, Mk0, EAxxac, ...
                EIxxac)
            % > Fern1und1991() finds SERRCalculator SERRCalculator for the J-integral and mode
            % partitioning only valid for a symmetric, balanced CLS joints,
            % where the adhesive interface is neglected in the lumped
            % overlap EA and EI
            %
            % > Source: [Fern1und1991]
            %------------------------------------------------------------
            
            % > Source: [Fern1und1991; EQ8]
            GI = Mk.^2/(4*EIxxac);
            
            % > Source: [Fern1und1991; EQ7]
            G_II_f = P.^2/(4*EAxxac);
            
            % > Source: [Fern1und1991; EQ9]
            G_II_b = (4*Mk.^2-Mk0.^2)./(16*EIxxac);
            
            % Add the Mode II contribution of the tensile and bending load
            GII = G_II_f+G_II_b;
            GT = GI+GII;
            MR = GII./GT;
        end
        
        function [GIdaf, GIIdaf, GTdaf, MRdaf] = SERRDAF(xserr, ...
                GI, GII, daf, x0daf)
            % > SERRDAF() re-calculates the given GI and GII on the
            % xserr locations under the influence of the DAF
            %
            %------------------------------------------------------------
            
            % The % SERR Footprint of the DAF is a set of x-coordinates
            % w.r.t. the DAF location. Using
            % linear interpolation, the % SERR Footprint of the DAF is
            % projected on the xserr coordinates of the intial SERR field
            %
            % Note: coordinates outside the x_daf will get 'NaN'
            
            % DAF % SERR footprint data x-coordinates
            xdaf = daf(:,1)+x0daf;
            
            % Interpolation of the DAF data points on the xserr coordinates
            dGI = interp1(xdaf, daf(:,2), xserr, 'linear');
            dGII = interp1(xdaf, daf(:,3), xserr, 'linear');
            
            % NaN values are outside the DAF data set; set to 0
            dGI(isnan(dGI))     = 0;
            dGII(isnan(dGII))   = 0;
            
            % Include the DAF effect
            GIdaf = GI+GI.*dGI;
            GIIdaf = GII+GII.*dGII;
            GTdaf = GIdaf+GIIdaf;
            MRdaf = GIIdaf./GTdaf;
        end
        
    end
end