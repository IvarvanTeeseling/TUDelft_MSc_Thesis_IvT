classdef clsOverlapLoads
    %clsOverlapLoads Summary of this class goes here
    %   Detailed explanation goes here
    
    % >> Source Bibliography
    % >
    % @book{daSilva2008,
    %   title={Modeling of adhesively bonded joints},
    %   author={Da Silva, Lucas Filipe Martins and {\"O}chsner, Andreas},
    %   year={2008},
    %   publisher={Springer}
    % }
    % >
    
    properties (GetAccess = public, SetAccess = private)
        Sxya    % Adhesive shear stress
        Syya    % Adhesive peel stress
        Nt      % N distribution - Top adherent (per unit width)
        Qt      % Q distribution - Top adherent (per unit width)
        Mt      % M distribution - Top adherent (per unit width)
        Nb      % N distribution - Bottom adherent (per unit width)
        Qb      % Q distribution - Bottom adherent (per unit width)
        Mb      % M distribution - Bottom adherent (per unit width)
    end
    
    methods
        function obj = clsOverlapLoads(x, Fk, Mk, Qk, lcb, Exx, ...
                t, Ea, Ga, ta)
            % > clsOverlapLoads() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before being stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Expected inputs for methodEng and sstate
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validnum        = @(x) isnumeric(x);
            
            % Add required arguments and their validation functions
            addRequired(p, 'x', validnum);
            addRequired(p, 'Fk', validnum);
            addRequired(p, 'Mk', validnum);
            addRequired(p, 'Qk', validnum);
            addRequired(p, 'lcb', validnum);
            addRequired(p, 'Exx', validnum);
            addRequired(p, 't', validnum);
            addRequired(p, 'Ea', validnum);
            addRequired(p, 'Ga', validnum);
            addRequired(p, 'ta', validnum);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            
            % Validate the input arguments
            parse(p, x, Fk, Mk, Qk, lcb, Exx, t, Ea, Ga, ta);
            
            % Adhesive Streses (from: Luo & Tong [2004, 2007])
            [Sxy, Syy] = clsOverlapLoads.adhesiveStress( ...
                p.Results.x, ...
                p.Results.Fk, ...
                p.Results.Mk, ...
                p.Results.Qk, ...
                p.Results.lcb, ...
                p.Results.Exx, ...
                p.Results.t, ...
                p.Results.Ea, ...
                p.Results.Ga, ...
                p.Results.ta);
            
            % Overlap adherent load distributions
            [Nt, Qt, Mt, Nb, Qb, Mb] = clsOverlapLoads.adherentLoads( ...
                p.Results.x(1,2)-p.Results.x(1,1), ...
                Sxy, ...
                Syy, ...
                p.Results.Fk, ...
                p.Results.Mk, ...
                p.Results.Qk, ...
                p.Results.t, ...
                p.Results.ta);
            
            % Perform convergence
            
            % Set object properties
            obj.Sxya = Sxy;
            obj.Syya = Syy;
            obj.Nt = Nt;
            obj.Qt = Qt;
            obj.Mt = Mt;
            obj.Nb = Nb;
            obj.Qb = Qb;
            obj.Mb = Mb;
        end
    end
    
    methods (Static)
        function [Sxy, Syy] = adhesiveStress(x, Fk, Mk, Qk, ...
                lcb, Exx, t, Ea, Ga, ta)
            % > adhesiveStress() find the adhesive peel (Syy) and
            % shear (Sxy) stress in the bonded overlap by applying the
            % Goland and Reissner (1944) solution framework
            %
            % Note: The equations are retreived from Luo and Tong,
            % published in 'Modeling of Adhesively Bonded Joints' as they
            % are more complete and do not contain the 'long joint'
            % simpliciation from Goland and Reissner in their bending
            % moment factor definition
            %
            %-------------------------------------------------------------
            
            % Peel stress integration constants
            %   > Source: [daSilva2008, page 33, EQ 2.20]
            beta_s = sqrt(2)/2*(24*Ea/(Exx*t^3*ta))^(1/4);
            beta_t = sqrt(8*Ga/(Exx*t*ta));
            
            % Integration constants
            %   > Source: [daSilva2008, page 42, EQ 2.65]
            Bs1 = Mk.*(sinh(beta_s*lcb).*cos(beta_s*lcb)+ ...
                cosh(beta_s*lcb).*sin(beta_s*lcb))+ ...
                Qk/beta_s.*sinh(beta_s*lcb).*sin(beta_s*lcb);
            Bs4 = Mk.*(sinh(beta_s*lcb).*cos(beta_s*lcb)- ...
                cosh(beta_s*lcb).*sin(beta_s*lcb))+ ...
                Qk/beta_s.*cosh(beta_s*lcb).*cos(beta_s*lcb);
            
            % Adhesive shear stress
            %   > Source: [daSilva2008, page 41, EQ 2.64]
            Sxy = beta_t*(Fk*t+6*Mk).*cosh(beta_t*x)./ ...
                (8*t*sinh(beta_t*lcb))+3*(Fk*t-2*Mk)./(8*t*lcb);
            
            % Adhesive peel stress
            %   > Source: [daSilva2008, page 41, EQ 2.64]
            Syy = 2*beta_s^2*(Bs1.*sinh(beta_s*x).*sin(beta_s*x)+ ...
                Bs4.*cosh(beta_s*x).*cos(beta_s*x))./ ...
                (sinh(2*beta_s*lcb)+sin(2*beta_s*lcb));
        end
        
        function [Nt, Qt, Mt, Nb, Qb, Mb] = adherentLoads(dx, ...
                Sxy, Syy, Fk, Mk, Qk, t, ta)
            % > adherentLoads() find the overlap (top and bottem)
            % adherent load distributions based on the force equilibrium of
            % an infinitesimal small element
            %
            % Note: integrations are performed numberical to avoid length,
            %       explicit equation implementation
            %
            % Note: It is important to use a sufficient amount of elements
            %       for an accurate represenation
            %
            %-------------------------------------------------------------
            
            % Top adherent - Axial force
            Nt = Fk-cumsum(Sxy.*dx,2,'omitnan');
            Nt(isnan(Sxy)) = NaN;
            
            % Top adherent - Shear force
            Qt = Qk-cumsum(Syy.*dx,2,'omitnan');
            Qt(isnan(Syy)) = NaN;
            
            % Top adherent - Bending moment
            Mt = Mk+cumsum(Qt*dx,2,'omitnan')-(t+ta)/2* ...
                cumsum(Sxy.*dx,2,'omitnan');
            Mt(isnan(Sxy)) = NaN;
            
            % Bottom adherent - Axial force
            Nb = cumsum(Sxy.*dx,2,'omitnan');
            Nb(isnan(Sxy)) = NaN;
            
            % Bottom adherent - Shear force
            Qb = cumsum(Syy.*dx,2,'omitnan');
            Qb(isnan(Syy)) = NaN;
            
            % Bottom adherent - Bending moment
            Mb = cumsum(Qb*dx,2,'omitnan')-(t+ta)/2* ...
                cumsum(Sxy.*dx,2,'omitnan');
            Mt(isnan(Sxy)) = NaN;
        end
        
    end
end

