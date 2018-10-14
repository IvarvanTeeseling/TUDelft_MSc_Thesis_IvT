classdef StressStrainConvert
    %StressStrainConvert Summary of this class goes here
    %   Detailed explanation goes here
    
    properties (GetAccess = public, SetAccess = private)
        
    end
   
    methods (Static)
        function econv = eConvert(e, direction)
            % > eConvert() converts engineering strain to true
            % strain and the other way around
            %------------------------------------------------------------
            
            % Input check
            if ~any(strcmp(direction, {'Eng2True' 'True2Eng'}))
                error('Wrong method input. Must be: Eng2True or True2Eng')
            end
            
            switch direction
                case 'Eng2True'
                    % Covert engineering strain to true strain
                    econv = log(1+e);
                case 'True2Eng'
                    % Convert true strain to engineering strain
                    econv = exp(e)-1;
            end
        end
        
        function sconv = sConvert(s, e, direction)
            % > sConvert() converts engineering stress to true
            % stress and the other way around
            %------------------------------------------------------------
            
            % Input check
            if ~any(strcmp(direction, {'Eng2True' 'True2Eng'}))
                error('Wrong method input. Must be: Eng2True or True2Eng')
            end
            
            switch direction
                case 'Eng2True'
                    % Covert engineering strain to true strain
                    sconv = s*(1+e);
                case 'True2Eng'
                    % Convert true strain to engineering strain
                    sconv = s./(1+e);
            end
        end
    end
end