classdef FRP < Material
    %FRP Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        
    end
    
    methods
        function obj = FRP(matLin, E1, E2, G, v12, varargin)
            % > FRP() constructs an instance of this class.
            %
            % > Inputs arguments are passed through an input parser to 
            % perform validity checks before begin stored as the instance
            % properties
            %-------------------------------------------------------------
            
            % Initialize an input parser scheme (input manager)
            p = inputParser;
            
            % Input argument validation functions
            validInpNum = @(x) isnumeric(x);
            validName   = @(x) ischar(x);
            
            % Add required arguments and their validation functions
            
            % Add optional arguments and their validation functions 
            % Note: optional argument default value = []
            addOptional(p, 'ct1', [], validInpNum);
            addOptional(p, 'ct2', [], validInpNum);
            addOptional(p, 'name', [], validName);
            
            % Validate the input arguments and store
            parse(p, varargin{:});  
            
            % Initialize inherited properties
            obj@Material(matLin, E1, E2, G, v12, ... 
                'ct1', p.Results.ct1, 'ct2', p.Results.ct2, ...
                'name', p.Results.name);
            
            % Set object properties
            obj.name = p.Results.name;
        end
    end
end

