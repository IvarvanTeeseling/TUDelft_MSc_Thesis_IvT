classdef MechStrainCycle
    %MechStrainCycle Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        exx0        % Laminate strain_xx cycle at z = 0 (x,y reference frame)
        kxx         % Laminate curvature_xx cycle (x,y reference frame)
        exx      	% Ply strain_xx cycle (x,y reference frame)
    end
    
    methods
        function obj = MechStrainCycle(zply, EAxx, EIxx, Nxx, Mxx, plyIdx)
            % > MechStrainCycle() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum = @(x) isnumeric(x);
            
            % Add required arguments and their validation functions
            addRequired(p, 'zply', validnum);
            addRequired(p, 'EAxx', validnum);
            addRequired(p, 'EIxx', validnum);
            addRequired(p, 'Nxx', validnum);
            addRequired(p, 'Mxx', validnum);
            addRequired(p, 'plyIdx', validnum);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            
            % Validate the input arguments
            parse(p, zply, EAxx, EIxx, Nxx, Mxx, plyIdx);
            
            % Input check spanning more than one paramater
            if any(plyIdx < 1) || plyIdx(end) > numel(zply)-1
                error('Ply index smaller than 1 or larger than the number of existing plies is not possible')
            end
            
            % Laminate membrane strain_xx and curvature_xx
            [exx0, kxx] = MechStrainCycle.laminateStrainCycle(Nxx, Mxx, ...
                EAxx, EIxx);
            
            % Ply strain_xx cycle for the requested plies
            exx = MechStrainCycle.plyStrainCycle(exx0, kxx, zply, plyIdx);
            
            % Set object properties
            obj.exx0 = exx0;
            obj.kxx = kxx;
            obj.exx = exx;
        end
    end
    
    methods (Static)
        function [ex0, kx] = laminateStrainCycle(Nxx, Mxx, EAxx, EIxx)
            % > laminateStrainCycle() find the laminate strain and curvature in
            % the xx direction as result from applied membrane load Nxx and
            % bending moment Mxx. Strain and curvature are calculated using
            % the plate smeared stiffness properties
            %
            % Note: smeared plate stiffness properties used as it is easier
            % to include the plane strain (through the width direction)
            %------------------------------------------------------------
            
            % Membrane strain and curvature
            ex0 = Nxx/EAxx;
            kx = Mxx/EIxx;
        end
        
        function exx = plyStrainCycle(ex0, kx, zply, plyIdx)
            % > laminateStrainCycle() find the ply strain cyle the xx
            % direction as result from the laminate strain and curvature
            % cycle
            %------------------------------------------------------------
            
            % Pre-allocate memory
            exx = zeros(size(kx));
            
            % Take strain from the central coordinate
            zmid = zply(1:end-1)-(zply(1:end-1)-zply(2:end))/2;
            
            % Initialize counter
            cnt = 1;
            
            for i = plyIdx
                % Note: exx is a [n m p q] matrix where:
                %   n = number of cracked elements + 1
                %   m = number of elements
                %   p = 2 (min and max load)
                %   q = number of plies
                
                % Total strain
                exx(:,:,:,cnt) = ex0-zmid(i)*kx;
                
                % Update counter
                cnt = cnt + 1;
            end
        end
    end  
end