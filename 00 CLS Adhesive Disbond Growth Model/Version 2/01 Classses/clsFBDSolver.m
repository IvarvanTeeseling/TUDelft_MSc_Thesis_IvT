classdef clsFBDSolver
    %clsFBDSolver Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        wac     % Vertical displacement (y-direction) in the AC section
        wcb     % Vertical displacement (y-direction) in the CB section
        Qac     % Shear distribution in the AC section
        Qcb     % Shear distribution in the CB section
        Mac     % Moment distribution in the AC section
        Mcb     % Moment distribution in the CB section
        Mk      % Mk = Mac(xac = lac0) ; Goland & Reissner
        Mk0     % Mk0 = Mcb(xcb = 0) ; Goland & Reissner
        Qk      % Qk = Qac(xac = lac0) ; Goland & Reissner
        Qk0     % Qk0 = Qcb(xcb = 0) ; Goland & Reissner
        Vk      % Qk based on the undeformed overlap
    end
    
    methods
        function obj = clsFBDSolver(xac, xcb, P, EIac, EIcb, ...
                lac, lcb, t, ta, BC)
            % > clsFBDSolver() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            expectedBC = {'RR' 'CC' 'RC'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum        = @(x) isnumeric(x);
            validBC         = @(x) any(strcmp(x, expectedBC));
            
            % Add required arguments and their validation functions
            addRequired(p, 'xac', validnum);
            addRequired(p, 'xcb', validnum);
            addRequired(p, 'P', validnum);
            addRequired(p, 'EIac', validnum);
            addRequired(p, 'EIcb', validnum);
            addRequired(p, 'lac', validnum);
            addRequired(p, 'lcb', validnum);
            addRequired(p, 't', validnum);
            addRequired(p, 'ta', validnum);
            addRequired(p, 'BC', validBC);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            
            % Validate the input arguments
            parse(p, xac, xcb, P, EIac, EIcb, lac, lcb, t, ta, BC);
            
            % Solve the CLS Free Body Diagram and find the overlap edge
            % loads by applying the Goland and Reissner Solution Framework
            GR = clsFBDSolver.overlapEdgeLoadsGR(p.Results.xac, ...
                p.Results.xcb, ...
                p.Results.P, ...
                p.Results.EIac, ...
                p.Results.EIcb, ...
                p.Results.lac, ...
                p.Results.lcb, ...
                p.Results.t, ...
                p.Results.ta, ...
                p.Results.BC);
            
            % Set object properties
            obj.wac = GR.wac;
            obj.wcb = GR.wcb;
            obj.Qac = GR.Qac;
            obj.Qcb = GR.Qcb;
            obj.Mac = GR.Mac;
            obj.Mcb = GR.Mcb;
            obj.Mk = GR.Mk;
            obj.Mk0 = GR.Mk0;
            obj.Qk = GR.Qk;
            obj.Qk0 = GR.Qk0;
            obj.Vk = GR.Vk;
        end
    end
    
    methods (Static)
        function GR = overlapEdgeLoadsGR(xac, xcb, P, EIac, EIcb, lac, ...
                lcb, t, ta, BC)
            % Solution eigenvalues
            lamba_ac = sqrt(P/EIac);
            lamba_cb = sqrt(P/EIcb);
            
            % Angle neutA1l axis w.r.t. adherent neutA1l axis
            alpha = (t+ta)./(2*(lac(1)+lcb(1)));
            
            % Find the diff. eq. integration constants
            switch BC
                case 'RR'
                    % Roller-Roller support boundary conditions
                    
                    % Reaction forces in A according to force equilibrium
                    RA = 0;
                    MA = 0;
                    
                    % Integration constants
                    A1 = 0;
                    B1 = -(sqrt(EIac)*(cosh(lamba_cb.*lcb)*t ...
                        +cosh(lamba_cb.*lcb)*ta+ ...
                        2*alpha*(lac+lcb)-t-ta))...
                        ./(2*(sqrt(EIac)*cosh(lamba_cb.*lcb).* ...
                        sinh(lamba_ac.*lac)+sqrt(EIcb).* ...
                        cosh(lamba_ac.*lac).*sinh(lamba_cb.*lcb)));
                    A0 = (-2*sqrt(EIcb/EIac)*B1.*cosh(lamba_ac.*lac).* ...
                        sinh(lamba_cb.*lcb)-2*alpha*(lac+lcb)+t+ta)./ ...
                        (2*cosh(lamba_cb.*lcb));
                    B0 = sqrt(EIcb/EIac)*B1.*cosh(lamba_ac.*lac);
                    
                case 'CC'
                    % Clamped-Clamped support boundary condition
                    % Note: point of applied load has a freedom of movement
                    % in the x-direction
                    
                    %Pre-allocate memory
                    RA = zeros(size(lcb,1), 1, size(P,2));
                    RB = zeros(size(lcb,1), 1, size(P,2));
                    MA = zeros(size(lcb,1), 1, size(P,2));
                    MB = zeros(size(lcb,1), 1, size(P,2));
                    A1 = zeros(size(lcb,1), 1, size(P,2));
                    B1 = zeros(size(lcb,1), 1, size(P,2));
                    A0 = zeros(size(lcb,1), 1, size(P,2));
                    B0 = zeros(size(lcb,1), 1, size(P,2));
                    
                    % Solve the linear set of 8 equations to find the Diff.
                    % Eq. integratino constants
                    for i = 1:size(P,3)
                        
                        % Applied load
                        Ptmp        = P(1,1,i);
                        lamba_actmp = lamba_ac(1,1,i);
                        lamba_cbtmp = lamba_cb(1,1,i);
                        
                        for j = 1:size(lcb,1)
                            % Free adherent (lac) and overlap (lcb) length
                            lcbtmp = lcb(j);
                            lactmp = lac(j);
                            
                            % Solve the linear system
                            amat = [0 0 0.1e1 / Ptmp 0 1 0 0 0; 0.1e1 / Ptmp 0 0 0 0 lamba_actmp 0 0; 0 0.1e1 / Ptmp * (lcbtmp + lactmp) 0.1e1 / Ptmp 0 0 0 cosh(lamba_cbtmp * lcbtmp) sinh(lamba_cbtmp * lcbtmp); 0 0.1e1 / Ptmp 0 0 0 0 lamba_cbtmp * sinh(lamba_cbtmp * lcbtmp) lamba_cbtmp * cosh(lamba_cbtmp * lcbtmp); 0.1e1 / Ptmp * lactmp -0.1e1 / Ptmp * lactmp 0 0 cosh(lamba_actmp * lactmp) sinh(lamba_actmp * lactmp) -1 0; 0.1e1 / Ptmp -0.1e1 / Ptmp 0 0 lamba_actmp * sinh(lamba_actmp * lactmp) lamba_actmp * cosh(lamba_actmp * lactmp) 0 -lamba_cbtmp; 1 -1 0 0 0 0 0 0; -lcbtmp - lactmp 0 -1 1 0 0 0 0];
                            dmat = [0 -alpha -(lcbtmp + lactmp) * alpha + t / 0.2e1 + ta / 0.2e1 -alpha -t / 0.2e1 - ta / 0.2e1 0 0 0]';
                            [L,U]=lu(amat);
                            bmat = U\(L\dmat);
                            
                            % Integration constants
                            RA(j,1,i) = bmat(1);
                            RB(j,1,i) = bmat(2);
                            MA(j,1,i) = bmat(3);
                            MB(j,1,i) = bmat(4);
                            A1(j,1,i) = bmat(5);
                            B1(j,1,i) = bmat(6);
                            A0(j,1,i) = bmat(7);
                            B0(j,1,i) = bmat(8);
                        end
                    end
            end
            
            % Overlap Edge Bending Moment
            GR.Mk = P.*(-A1.*cosh(lamba_ac.*lac)- ...
                B1.*sinh(lamba_ac.*lac));
            % Overlap Edge Shear Force
            GR.Qk = P.*(-A1.*lamba_ac.*sinh(lamba_ac.*lac)- ...
                B1.*lamba_ac.*cosh(lamba_ac.*lac));
            % Overlap Lumped Edge Bending Moment
            GR.Mk0 = P.*(-A0.*cosh(lamba_cb.*0)- ...
                B0.*sinh(lamba_cb.*0));
            % Overlap Lumped Edge Bending Moment
            GR.Qk0 = P.*(-A0.*lamba_cb.*sinh(lamba_cb.*0)- ...
                B0.*lamba_cb.*cosh(lamba_cb.*0));
            
            % Vk = Qk according to force equilibrium of the undeformed body
            % Note: Qk is found based on the geometric non-linear solution
            %       of the Free Body Diagram; it derrived from the force
            %       equilibirum after significant deformation
            %
            %       To solve the adhesive stresses, Mk and Qk will be
            %       applied as the overlap edge loads. However, during this
            %       step, the overlap is assumed to be geometric
            %       non-linear. In other words: the deformed state on which
            %       Qk has been computed does not exist anymore. To account
            %       for this change many solutions adopt Vk in their
            %       adhesive stress formulation instead of Qk.
            GR.Vk = (P*(t+ta)-2*GR.Mk)./(2*lcb);
            
            % Vertical displacement - Free adherent section
            GR.wac = A1.*cosh(lamba_ac.*xac)+ ...
                B1.*sinh(lamba_ac.*xac)+(alpha+RA./P).*xac+MA./P;
            % Vertical displacement - Lumped overlap section
            GR.wcb = A0.*cosh(lamba_cb.*xcb)+ ...
                B0.*sinh(lamba_cb.*xcb)+alpha.*(lac+xcb)- ...
                (t+ta)/2+(lac+xcb).*RA./P+MA./P;
            
            % Moment distribution - Free adherent section
            GR.Mac = P.*(-A1.*cosh(lamba_ac.*xac)- ...
                B1.*sinh(lamba_ac.*xac));
            % Moment distribution - Lumped overlap section
            GR.Qac = P.*(-A1.*lamba_ac.*sinh(lamba_ac.*xac)- ...
                B1.*lamba_ac.*cosh(lamba_ac.*xac));
            
            % Shear distribution - Free adherent section
            GR.Mcb = P.*(-A0.*cosh(lamba_cb.*xcb)- ...
                B0.*sinh(lamba_cb.*xcb));
            % Shear distribution - Lumped overlap section
            GR.Qcb = P.*(-A0.*lamba_cb.*sinh(lamba_cb.*xcb)- ...
                B0.*lamba_cb.*cosh(lamba_cb.*xcb));
        end
    end
end

