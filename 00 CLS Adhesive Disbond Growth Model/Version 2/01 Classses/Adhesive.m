classdef Adhesive < Material
    %Adhesive Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        c0      % D. Burger Paris-Law C coefficient at MR = 0
        c100    % D. Burger Paris-Law C coefficient at MR = 100
        m0      % D. Burger Paris-Law m coefficient at MR = 0
    end
    
    methods
        function obj = Adhesive(matLin, E1, E2, G, v12, c0, ... 
                c100, m0, varargin)
            % > Adhesive() constructs an instance of this class.
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
            addRequired(p, 'c0', validInpNum);
            addRequired(p, 'c100', validInpNum);
            addRequired(p, 'm0', validInpNum);
            addOptional(p, 'ct1', [], validInpNum);
            addOptional(p, 'ct2', [], validInpNum);
            addOptional(p, 'name', [], validName);
            
            % Validate the input arguments
            parse(p, c0, c100, m0, varargin{:});
            
            % Initialize inherited properties
            obj@Material(matLin, E1, E2, G, v12, ...
                'ct1', p.Results.ct1, 'ct2', p.Results.ct2, ...
                'name', p.Results.name);
            
            % Set object properties
            obj.c0 = p.Results.c0;
            obj.c100 = p.Results.c100;
            obj.m0 = p.Results.m0;
            obj.name = p.Results.name;
        end
    end
end

