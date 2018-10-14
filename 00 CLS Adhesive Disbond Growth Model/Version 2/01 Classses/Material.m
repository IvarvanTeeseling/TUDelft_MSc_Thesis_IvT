classdef Material
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        matLin      % Material linearity
        E1          % Young's Modulus 1-direction
        E2          % Young's Modulus 2-direction
        G           % Shear modulus
        v12         % Poisson ratio 12-direction
        ct1         % Thermal expansion coeficient 1-direction (optional)
        ct2         % Thermal expansion coeficient 2-direction (optional)
        name        % Material name (optional)
    end
    
    methods
        function obj = Material(matLin, E1, E2, G, v12, varargin)
            % > Material() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to 
            % perform validity checks before begin stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % matLin must either be 'Elastic' or 'Plastic' 
            expectedmatLin = {'elastic' 'plastic'};

            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validInpVal = @(x) isnumeric(x);
            validmatLin = @(x) any(strcmp(x, expectedmatLin)); 
            validName   = @(x) ischar(x);
            
            % Add required arguments and their validation functions
            addRequired(p, 'matLin', validmatLin)
            addRequired(p, 'E1', validInpVal);
            addRequired(p, 'E2', validInpVal);
            addRequired(p, 'G', validInpVal);
            addRequired(p, 'v12', validInpVal);
            
            % Add optional arguments and their validation functions 
            % Note: optional argument default value = []
            addOptional(p, 'ct1', [], validInpVal);
            addOptional(p, 'ct2', [], validInpVal);
            addOptional(p, 'name', [], validName);
            
            % Validate the input arguments
            parse(p, matLin, E1, E2, G, v12, varargin{:});  
            
            % Set object properties
            obj.matLin = p.Results.matLin;
            obj.E1 = p.Results.E1;
            obj.E2 = p.Results.E2;
            obj.G = p.Results.G;
            obj.v12 = p.Results.v12;
            obj.ct1 = p.Results.ct1;
            obj.ct2 = p.Results.ct2;
            obj.name = p.Results.name;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

