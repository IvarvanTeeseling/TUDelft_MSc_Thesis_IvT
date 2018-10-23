classdef MetalStrainCycle
    %MetalStrainCycle Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        exx         % Metal ply total (true/eng) strain xx
        exxmech   	% Metal ply mechanical (true/eng) strain xx
        exxpl     	% Metal ply plastic (true/eng) strain xx
        exxres   	% Metal ply residual (true/eng) strain xx
        eTypeInp    % Input type (default = 'EngStrain') (optional)
        eTypeOut    % Output type (default = 'EngStrain') (optional)
    end
    
    methods
        function obj = MetalStrainCycle(Exx, exxmech, exxres, exxrth, varargin)
            % > MetalStrainCycle() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            expectedeMethod = {'CutOff' 'TrueStressStrain'};
            expectedeType = {'TrueStrain' 'EngStrain'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum       	= @(x) isnumeric(x);
            valideMethod    = @(x) any(strcmp(x, expectedeMethod));
            valideType      = @(x) any(strcmp(x, expectedeType));
            
            % Add required arguments and their validation functions
            addRequired(p, 'Exx', validnum);
            addRequired(p, 'exxm', validnum);
            addRequired(p, 'exxr', validnum);
            addRequired(p, 'exxrth', validnum);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            %       1. eMethod = noeMethod
            addOptional(p, 'eMethod', 'Material', valideMethod)
            addOptional(p, 'eTypeInp', 'EngStrain', valideType)
            addOptional(p, 'eTypeOut', 'EngStrain', valideType)
            
            % Validate the input arguments
            parse(p, Exx, exxmech, exxres, exxrth, varargin{:});
            
            obj.eTypeInp = p.Results.eTypeInp;
            obj.eTypeOut = p.Results.eTypeOut;
            
            % Input check spanning more than one paramater
            
            if strcmp(p.Results.eTypeInp, 'EngStrain')
                % Input is in engineering strain > convert
                exxm_tr = StressStrainConvert.eConvert(exxmech, 'Eng2True');
                exxr_tr = StressStrainConvert.eConvert(exxres, 'Eng2True');
                exxrth_tr = StressStrainConvert.eConvert(exxrth, 'Eng2True');
            else
                % Input is already true strain
                exxm_tr = exxmech;
                exxr_tr = exxres;
                exxrth_tr = exxrth;
            end
            
            % Elastic Young's Modulus from material data
            Eel = (Exx(2,2)-Exx(1,2))/(Exx(2,1)-Exx(1,1));
            
            % Flag for residual strain (1 = On, 0 = Off)
            flag = ones(size(exxm_tr));
            
            % Pre-allocate
            exx_tr = zeros(size(exxm_tr));
            exxpl_tr = zeros(size(exxm_tr));
            
            switch p.Results.eMethod
                case 'CutOff'
                    % Find the total strain
                    exx_tr = exxm_tr+exxr_tr(1);
                    
                    % Strain cut off exceeded the yield strain
                    idx = exx_tr > (Exx(2,1)-Exx(1,1));
                    
                    % Assumed: cut-off strain equals plastic strain
                    exxpl_tr(idx) = exx_tr(idx)-(Exx(2,1)-Exx(1,1));
                    
                    % Residual strain is reversed during plasticity
                    exx_tr(idx) = (Exx(2,1)-Exx(1,1));
                    
                case 'TrueStressStrain'
                    for i = 1:size(exxm_tr, 1)
                        % Total element true strain
                        exx_tr(i,:,:,:) = exxm_tr(i,:,:,:) + ...
                            exxr_tr(1)*flag(i,:,:,:);
                        
                        % Total element true stress
                        Sxxtot_tr = MetalStressCycle.findTrueStress( ...
                            exx_tr(i,:,:,:), Exx);
                        
                        % Find the element plastic strain
                        epl_tr_new = exx_tr(i,:,:,:)- ...
                            Sxxtot_tr/Eel;
                        
                        % Remove floating point errors
                        epl_tr_new(abs(epl_tr_new)<1e-8) = 0;

                        exxpl_tr(i,:,:,:) = epl_tr_new;
                        
                        % Removed reversed residual strain
                        idx = exxpl_tr(i,:,:,:) > exxrth_tr;
                        flag(i:end, idx==1) = 0;
                        
                        % Re-assemble the final total true strain cycle
                        exx_tr(i,:,:,:) = exxm_tr(i,:,:,:) + ...
                            exxr_tr(1)*flag(i,:,:,:);
                    end
            end
            
            % Set output to the requested format
            if strcmp(p.Results.eTypeOut, 'EngStrain')
                obj.exx = StressStrainConvert.eConvert(exx_tr, ... 
                    'True2Eng');
                obj.exxmech = StressStrainConvert.eConvert(exxm_tr, ...
                    'True2Eng');
                obj.exxpl = StressStrainConvert.eConvert(exxpl_tr, ...
                    'True2Eng');
                obj.exxres =StressStrainConvert.eConvert(exxr_tr(1)*flag, ...
                    'True2Eng');
            else
                obj.exx = exx_tr;
                obj.exxmech = exxm_tr;
                obj.exxpl = exxpl_tr;
                obj.exxres = exxr_tr(1)*flag;
            end
        end
    end
    
    methods (Static)
        function econv = eConvert(e, direction)
            % > eConvert() converts engineering strain to true
            % strain and the other way around
            %------------------------------------------------------------
            
            % Input check
            if ~any(strcmp(direction, {'Eng2True' 'True2Eng'}))
                error('Wrong method input. Must be: Eng2True or True2Eng')
            end
            
            switch direction
                case 'Eng2True'
                    % Covert engineering strain to true strain
                    econv = log(1+e);
                case 'True2Eng'
                    % Convert true strain to engineering strain
                    econv = exp(e)-1;
            end
        end
    end
end