function [x12] = xMatrix_Zipper(x1, x2, val, option)

if ~isnumeric(val)
    error('xMatrix_Zipper: <val> input must be numeric.')
end

if ~any(strcmp(option, {'zip' 'unzip'}))
    error('xMatrix_Zipper: invalid input for <option>')
end

switch option
    case 'zip'
        if isnan(val)
            % Non active elements indicated by NaN
            x1(isnan(x1)) = 0;
            x2(isnan(x2)) = 0;
        elseif isinf(val)
            % Non active elements indicated by Inf
            x1(isinf(x1)) = 0;
            x2(isinf(x2)) = 0;
        else
            % Non active elements indicated by a value
            x1(x1==val) = 0;
            x2(x1==val) = 0;
        end
        
        % Disbond increments
        qc = size(x1, 1)-1;
        
        % Zip together
        x12 = [x1(:,1:end-qc,:) x1(:,end-qc+1:end,:)+x2(:,1:qc,:) x2(:,qc+1:end,:)];
    case 'unzip'

end

end