classdef AdherentStiffness
    %AdherentStiffness Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        EAxxac
        EIxxac
        EAxxcb
        EIxxcb
    end
    
    methods
        function obj = AdherentStiffness(Exx, Exxa, t, ta)
            % > AdherentStiffness() constructs an instance of this class.
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
            addRequired(p, 'Exx', validnum);
            addRequired(p, 'Exxa', validnum);
            addRequired(p, 't', validnum);
            addRequired(p, 'ta', validnum);
            
            % Add optional arguments and their validation functions
            % Note: optional argument default values:
            
            % Validate the input arguments
            parse(p, Exx, Exxa, t, ta);
            
            % Get axial stiffness
            EAxxac = Exx*t;
            EAxxcb = Exx*t*2 + Exxa*ta;
            
            % Get bending stiffness
            EIxxac = Exx*1/12*t^3;
            EIxxcb = Exx*2*(1/12*t^3+t*(t/2+ta/2)^2) + Exxa*1/12*ta^3/12;
            
            % Set object properties
            obj.EAxxac = EAxxac;
            obj.EIxxac = EIxxac;
            obj.EAxxcb = EAxxcb;
            obj.EIxxcb = EIxxcb;
        end
    end
end

