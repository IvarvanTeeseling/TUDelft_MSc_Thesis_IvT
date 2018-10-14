classdef clsAdherent
    %clsAdherent Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        L           % Total length
        q           % Number of elements
        x           % Element x-coordinates
        t           % Thickness
        xpos        % Default = 'Central'
    end
    
    methods
        function obj = clsAdherent(L, t, q, varargin)
            % > clsAdherent() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            expectedxpos = {'LBoundary', 'Central', 'RBoundary'};
            expectedrfpos = {'Left' 'Right'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum        = @(x) isnumeric(lac0);
            validq          = @(x) isnumeric(q) && floor(x) == x;
            validxpos       = @(x) any(strcmp(x, expectedxpos));
            validrfpos      = @(x) any(strcmp(x, expectedrfpos));
            
            % Add required arguments and their validation functions
            addRequired(p, 'L', validnum);
            addRequired(p, 't', validnum);
            addRequired(p, 'Exx', validnum);
            addRequired(p, 'q', validq);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            %       1. qcrack = 0
            %       2. Method = 'Central'
            %       3. rfpos = 'Left'
            %       3. pasval = NaN
            addOptional(p, 'qcrack', 0, validqcrack);
            addOptional(p, 'method', 'Central', validxpos);
            addOptional(p, 'rfpos', 'Left', validrfpos);
            addOptional(p, 'pasval', NaN, validnum);
            
            % Validate the input arguments
            parse(p, lac0, lcb0, q, varargin{:});
            
            % Discretize the adherent and return element x-coordinates
            Xvec = discretizeAdherent(p.results.L, p.results.q, ...
                p.results.rfpos, p.results.xpos);
            
            % Add add crack increments (represented by the rows)
            Xmat = includeCrackElements(Xvec, p.results.qcrack, ...
                p.results.pasval);
            
            % Set object properties
            
        end
    end
    
    methods (Static)
        function xvec = discretizeAdherent(L, q, rfpos, xpos)
            % > discretizeAdherent() discretizes the adherent into q
            % elements and defines the x-coordinates of each element based
            % on the reference frame position and the x-position within an
            % element
            %
            % > Several functions are called that are defined in this
            % static method domain of this class
            %------------------------------------------------------------
            
            % Element length
            dL = L/q;
            
            % Reference frame (x, y) position
            switch rfpos
                case 'Left'
                    xleft = 0;
                case 'Right'
                    xleft = -L;
            end
            
            % Discretize and return x-coordinates depending on the x
            % location with each element
            switch xpos
                case 'LBoundary'
                    xvec = xleft:dL:(L+xleft)-dL;
                case 'Central'
                    xvec = xleft+dL/2:dL:(L+xleft)-dL/2;
                case 'RBoundary'
                    xvec = xleft+dL:dL:(L+xleft);
            end
        end
        
        function [xmat1, xmat2, Lvec1, Lvec2] = includeCrackElements(...
                xvec1, L1, xvec2, L2, qcrack, pasval)
            % > includeCrackElements() adds cracked increments to form a
            % matrix 'xmat' where each row represents the updated adherent
            % geometry and corresponding x-coordinates after an element has
            % been cracked
            %
            % Elements are removed from xvec2 and added to xvec1
            %
            % > Several functions are called that are defined in this
            % static method domain of this class
            %------------------------------------------------------------ 
            
            % Length of a crack-element
            % Note:
            %   > Elements from xvec2 are cracked, removed and added to
            %   xvec1
            dL2 = xvec2(2)-xvec2(1);
            
            % Store adherent length after each crack increment
            Lvec1 = L1:dL2:L1+qcrack*dL2;       % Becomes longer
            Lvec2 = L2:-dL2:L2-qcrack*dL2;      % Becomes shorter
            
            % Pre-fill with pasval which flags passive elements
            xmat1 = ones(qcrack+1, size(xvec1, 2)+qcrack)*pasval;
            xmat2 = ones(qcrack+1, size(xvec2, 2))*pasval;
            
            %
            xmat1(1,1:end-qrack) = xvec1;
            xmat2(1,1+qcrack:end) = xvec2;
            
            for i = 1:qcrack
                if any(xvec1 < 0)
                    % Reference frame on the right boundary
                    xmat1(i+1,1:size(xvec1+i)) = [xmat1(1     
                else
                    % Reference frame on the left boundary
                end
            end
            
        end
        
    end
end

