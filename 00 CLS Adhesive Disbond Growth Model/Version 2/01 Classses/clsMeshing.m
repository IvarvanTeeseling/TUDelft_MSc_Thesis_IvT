classdef clsMeshing
    %clsMeshing Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        xAC         % x-coordinates (AC reference frame)
        xCB         % x-coordinates (CB reference frame)
        xAB         % x-coordinates (AB reference frame)
        xBC         % x-coordinates (BC reference frame)
        lAC         % Adherent length section AC
        lCB         % Adherent length section CB
        b           % Disbond incremental length (optional)
        xLoc        % Location for the element x-coordinate (optional)
        pasVal      % Value for passive elements (optional)
    end
    
    methods
        function obj = clsMeshing(lac0, lcb0, q, varargin)
            % > clsMeshing() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            expectedxloc = {'L', 'C', 'R'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum        = @(x) isnumeric(x);
            validq          = @(x) all(isnumeric(x) & floor(x) == x & ... 
                (size(x, 2) == 1 || size(x, 2) == 2));
            validxloc       = @(x) any(strcmp(x, expectedxloc));
            
            % Add required arguments and their validation functions
            addRequired(p, 'lac0', validnum);
            addRequired(p, 'lcb0', validnum);
            addRequired(p, 'q', validq);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            %       1. qcrack = 0
            %       2. xloc = 'Central'
            %       3. pasval = NaN
            addOptional(p, 'qcrack', 0, validq);
            addOptional(p, 'xloc', 'Central', validxloc);
            addOptional(p, 'pasval', NaN, validnum);
            
            % Validate the input arguments
            parse(p, lac0, lcb0, q, varargin{:});
            
            % Repeat 'q' for AC and BC if only one 'q' is given
            if size(p.Results.q, 2) == 1
                p.Results.q = [p.Results.q p.Results.q];
            end
            
            % Section x-coordinates
            xac = clsMeshing.discretizeAdherent(p.Results.lac0, ...
                p.Results.q(1), 'Left', p.Results.xloc);
            xcb = clsMeshing.discretizeAdherent(p.Results.lcb0, ...
                p.Results.q(2), 'Left', p.Results.xloc);
            xbc = clsMeshing.discretizeAdherent(p.Results.lcb0, ...
                p.Results.q(2), 'Right', p.Results.xloc);
            xab = [xac xcb+lac0];
            
            % Expand with cracked elements
            [xAC, xCB, lAC, lCB, b] = clsMeshing.crackedElements(...
                xac, xcb, p.Results.lac0, p.Results.lcb0, ...
                p.Results.qcrack, p.Results.pasval);
            
            [~, xBC, ~, ~, ~] = clsMeshing.crackedElements(...
                xac, xbc, p.Results.lac0, p.Results.lcb0, ...
                p.Results.qcrack, p.Results.pasval);

            % Set object properties
            obj.xAC = xAC;
            obj.xCB = xCB;
            obj.xBC = xBC;
            obj.xAB = xab;
            obj.lAC = lAC;
            obj.lCB = lCB;
            obj.b = b;
            obj.xLoc = p.Results.xloc;
            obj.pasVal = p.Results.pasval;
        end
    end
    
    methods (Static)
        function xvec = discretizeAdherent(l, q, rfpos, xpos)
            % > discretizeAdherent() discretizes the adherent into q
            % elements and defines the x-coordinates of each element based
            % on the reference frame position and the x-position within an
            % element
            %
            % > Several functions are called that are defined in this
            % static method domain of this class
            %------------------------------------------------------------
            
            % Input checking
            if ~any(strcmp(rfpos, {'Left' 'Right'}))
                error('Invalid reference frame position.')
            end
            if ~any(strcmp(xpos, {'L' 'C' 'R'}))
                error('Invalid element x-position.')
            end
            
            % Element length
            dl = l/q;
            
            % Reference frame (x, y) position
            switch rfpos
                case 'Left'
                    xleft = 0;
                case 'Right'
                    xleft = -l;
            end
            
            % Discretize and return x-coordinates depending on location
            % within each element to set the x-coordinate
            switch xpos
                case 'L'
                    % x at taken from element left boundary location
                    xvec = xleft:dl:(l+xleft)-dl;
                case 'C'
                    % x at taken from element central location
                    xvec = xleft+dl/2:dl:(l+xleft)-dl/2;
                case 'R'
                    % x at taken from element right boundary location
                    xvec = xleft+dl:dl:(l+xleft);
            end
        end
        
        function [xAC, xCB, lac, lcb, b] = crackedElements(xac, xcb, ... 
                lac0, lcb0, qc, pv)
            % > discretizeAdherent()
            %
            % > Several functions are called that are defined in this
            % static method domain of this class
            %------------------------------------------------------------
            
            % Disbond increment element length (elements from cb region)
            dlcb = xcb(2)-xcb(1);
            
            % Disbond increment elements get added to section AC
            lac = (lac0:dlcb:lac0+qc*dlcb)';
            
            % Disbond increment elements get removed from section BC
            lcb = (lcb0:-dlcb:lcb0-dlcb*qc)';
            
            % Disbond increments
            b = dlcb*(0:1:qc);
            
            % Reference frame must be at the left hand edge for further
            % computations; adjust if located on the right hand edge
            if any (xac < 0)
                % AC reference frame on right edge; adjustment needed
                xac_adj = xac+lac0;
                xacFlag = 1;
            else
                % AC reference frame on the left edge; no adjustement
                xac_adj = xac;
                xacFlag = 0;
            end
            if any(xcb < 0)
                % CB reference frame on the right edge; adjustment needed
                xcb_adj = xcb+lcb0;
                xcbFlag = 1;
            else
                % CB reference frame on the left edge; no adjustement
                xcb_adj = xcb;
                xcbFlag = 0;
            end
            
            % Create x-vector spanning the entire CLS (from A to B)
            xab = [xac_adj(1,:) xcb_adj(1,:)+lac0];
            
            % Pre-allocate passive elements (elements that will be
            % added/removed during future disbond increments) with their
            % marker value for identification
            xAC = ones(1+qc, size(xac,2)+qc)*pv;
            xCB = ones(1+qc, size(xcb,2))*pv;
            
            % Add/remove elements for each disbond increment
            for i = 1:qc+1
                idx = size(xac,2)+(i-1);
                % dCB element added to section AC
                xAC(i, 1:idx) = xab(1, 1:idx)-xacFlag*(lac(i));
                % dCB element removed from section CB
                xCB(i, i:end) = xcb_adj(1:end-(i-1))-xcbFlag*lcb(i);
            end
        end
    end
end

