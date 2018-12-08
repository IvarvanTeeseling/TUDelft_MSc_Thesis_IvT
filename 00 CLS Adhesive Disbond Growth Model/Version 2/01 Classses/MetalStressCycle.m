classdef MetalStressCycle
    %MetalStressCycle Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        Sxx         % Metal ply total stress (True/Eng)
        Sxxm       % Metal ply cycle mean stress (True/Eng)
        Sxxa       % Metal ply cycle amplitude stress (True/Eng)
        R           % Metal ply load ratio (Smin/Smax) (True/Eng)
        eTypeIn 	% Input strain type (default = 'EngStrain') (optional)
        sTypeOut 	% Output stress type (default = 'EngStress') (optional)
    end
    
    methods
        function obj = MetalStressCycle(Exx, exx, varargin)
            % > MetalStressCycle() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            expectedeType = {'TrueStrain' 'EngStrain'};
            expectedsType = {'TrueStress' 'EngStress'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum        = @(x) isnumeric(x);
            valideTypeIn    = @(x) any(strcmp(x, expectedeType));
            validsType      = @(x) any(strcmp(x, expectedsType));
            
            % Add required arguments and their validation functions
            addRequired(p, 'E', validnum);
            addRequired(p, 'exx', validnum);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            %       1. eTypeIn = 'EngStrain'
            %       1. sTypeOut = 'EngStess'
            addOptional(p, 'eTypeIn', 'EngStrain', valideTypeIn);
            addOptional(p, 'sTypeOut', 'EngStess', validsType);
            
            % Validate the input arguments
            parse(p, Exx, exx, varargin{:});
            
            % Input check spanning more than one paramater
            
            obj.eTypeIn = p.Results.eTypeIn;
            obj.sTypeOut = p.Results.sTypeOut;
            
            if strcmp(p.Results.eTypeIn, 'EngStrain')
                % Input is Engineering Strain > convert to True Strain
                exx_tr = MetalStrainCycle.strainConverter(exx, 'Eng2True');
            else
                % Input is True Strain > no convert
                exx_tr = exx;
            end
            
            % Find the metal ply stresses using the stress-strain curve
            Sxx_tr = MetalStressCycle.findTrueStress(exx_tr, Exx);
            
            if strcmp(p.Results.sTypeOut, 'EngStress')
                % Output is Engineering Stress > convert
                Sxx = StressStrainConvert.sConvert(Sxx_tr, exx, 'True2Eng');
            else
                % Output is True Stress > no convert
                Sxx = Sxx_tr;
            end
            
            % Find the load ratio
            R = Sxx(:,:,1,:)./Sxx(:,:,2,:);
            
            % Find the cycle mean stress
            obj.Sxxm = (1+R)./2.*Sxx(:,:,2,:);
            
            % Find the cycle amplitude stress
            obj.Sxxa = (1-R)./2.*Sxx(:,:,2,:);
            
            % Set output
            obj.Sxx = Sxx;
            obj.R = R;
            
        end
    end
    
    methods (Static)
        function s_tr = findTrueStress(e_tr, E)
            % > findTrueStress() finds the true stress based on the
            % material elastic-plastic data given
            %------------------------------------------------------------
            
            % Input check
            if size(E, 2) ~= 2 || size(E, 1) < 3
                error('E must have 2 columns and a minimum of 3 rows')
            end
            
            % Use interpolation to find true stress
            % Note: data outside the E range will get 'NaN' as value
            s_tr = interp1(E(:,1), E(:,2), e_tr, 'linear');
        end
    end
end
