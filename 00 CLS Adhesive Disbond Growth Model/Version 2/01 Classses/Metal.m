classdef Metal < Material
    %Metal Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        sy1     % Yield strength 1-direction (optional)
        sy2     % Yield strength 2-direction (optional) 
        su1     % Ultimate strength 1-direction (optional)
        su2     % Ultimate strength 2-direction (optional)
    end
    
    methods
        function obj = Metal(matLin, E1, E2, G, v12, varargin)
            % > Metal() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to 
            % perform validity checks before begin stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validInpNum = @(x) isnumeric(x);
            validSuSy   = @(x) isnumeric(x) && (x > 0);
            validName   = @(x) ischar(x);
            
            % Add required arguments and their validation functions
            
            % Add optional arguments and their validation functions 
            % Note: optional argument default value = []
            addOptional(p, 'sy1', [], validSuSy); 
            addOptional(p, 'sy2', [], validSuSy);   
            addOptional(p, 'su1', [], validSuSy); 
            addOptional(p, 'su2', [], validSuSy);
            addOptional(p, 'ct1', [], validInpNum);
            addOptional(p, 'ct2', [], validInpNum);
            addOptional(p, 'name', [], validName);
            
            % Validate the input arguments
            parse(p, varargin{:});  
            
            % Initialize inherited properties
            obj@Material(matLin, E1, E2, G, v12, ... 
                'ct1', p.Results.ct1, 'ct2', p.Results.ct2, ...
                'name', p.Results.name);
            
            % Set object properties
            obj.sy1 = p.Results.sy1;
            obj.sy2 = p.Results.sy2;
            obj.su1 = p.Results.su1;
            obj.su2 = p.Results.su2;
            obj.name = p.Results.name;
        end
    end
end

