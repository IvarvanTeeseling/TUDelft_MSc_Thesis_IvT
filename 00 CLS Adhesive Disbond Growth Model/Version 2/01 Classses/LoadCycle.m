classdef LoadCycle
    % LoadCycle Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (Dependent)
        P           % Load cycle vector in [N]
        Prunning    % Load cycle vector in [N/m]
        S           % Load cycle vector in [N/m^2]
    end
    
    properties (GetAccess = private, Hidden)
        inpFType    % Load input type; [N], [N/m] or [N/m^2]
        inpF        % Maximum load
        inpR        % R-ratio
        inpd        % Specimen width
        inpt        % Specimen thickness
    end
    
    methods
        function obj = LoadCycle(FType, F, R, d, t)
            % > LoadCycle() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to
            % perform validity checks before begin stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % FType must either 'N', 'N/m' or 'N/m^2'
            expectedFType = {'N' 'N/m' 'N/m^2'};
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validFType  = @(x) any(strcmp(x, expectedFType));
            validInpNum = @(x) isnumeric(x);
            validR      = @(x) isnumeric(x) && (x<1) && (x>-1);
            
            % Add required arguments and their validation functions
            addRequired(p, 'FType', validFType)
            addRequired(p, 'F', validInpNum)
            addRequired(p, 'R', validR)
            addRequired(p, 'd', validInpNum)
            addRequired(p, 't', validInpNum)
            
            % Add optional arguments and their validation functions
            % Note: optional argument default value = []
            
            % Validate the input arguments
            parse(p, FType, F, R, d, t);
            
            obj.inpFType = p.Results.FType;
            obj.inpF = p.Results.F;
            obj.inpR = p.Results.R;
            obj.inpd = p.Results.d;
            obj.inpt = p.Results.t;
        end
        
        function P = get.P(obj)
            % > get.P(obj) performs a mathmetical operation to determine
            % the class property P when requested
            %-------------------------------------------------------------
            switch obj.inpFType
                case 'N'
                    P(1,1,1) = obj.inpF*obj.inpR;
                    P(1,1,2) = obj.inpF;
                case 'N/m'
                    P(1,1,1) = obj.inpF*obj.inpd*obj.inpR;
                    P(1,1,2) = obj.inpF*obj.inpd;
                case 'N/m^2'
                    P(1,1,1) = obj.inpF*obj.inpd*obj.inpt*obj.inpR;
                    P(1,1,2) = obj.inpF*obj.inpd*obj.inpt;
            end
        end
        
        function Prunning = get.Prunning(obj)
            % > get.Prunning(obj) performs a mathmetical operation to
            % determine the class property Prunning when requested
            %-------------------------------------------------------------
            switch obj.inpFType
                case 'N'
                    Prunning(1,1,1) = obj.inpF/obj.inpd*obj.inpR;
                    Prunning(1,1,2) = obj.inpF/obj.inpd;
                case 'N/m'
                    Prunning(1,1,1) = obj.inpF*obj.inpR;
                    Prunning(1,1,2) = obj.inpF;
                case 'N/m^2'
                    Prunning(1,1,1) = obj.inpF*obj.inpt*obj.inpR;
                    Prunning(1,1,2) = obj.inpF*obj.inpt;
            end
        end
        
        function S = get.S(obj)
            % > get.S(obj) performs a mathmetical operation to determine
            % the class property S when requested
            %-------------------------------------------------------------
            switch obj.inpFType
                case 'N'
                    S(:,:,1) = obj.inpF/(obj.inpd*obj.inpt)*obj.inpR;
                    S(:,:,2) = obj.inpF/(obj.inpd*obj.inpt);
                case 'N/m'
                    S(:,:,1) = obj.inpF/obj.inpd*obj.inpR;
                    S(:,:,2) = obj.inpF/obj.inpd;
                case 'N/m^2'
                    S(:,:,1) = obj.inpF*obj.inpR;
                    S(:,:,2) = obj.inpF;
            end
        end
    end
end

