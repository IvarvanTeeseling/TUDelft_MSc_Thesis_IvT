function [xA, xAB, xB, xBB, l_A, l_B] = Discretize_Adherents(l_A0, l_B0, q, q_max, val, method)

if ~isnumeric(val)
    error('DISCRETIZE_ADHERENTS: val must be numeric.')
end

if ~any(strcmp(method, {'LeftBoundary', 'Central', 'RightBoundary'}))
    error('DISCRETIZE_ADHERENTS: incorrect method entered.')
end
    
% Element length
dlA = l_A0/q;
dlB = l_B0/q;

switch method
    case 'LeftBoundary'
        xAtmp = 0:dlA:l_A0-dlA;
        xBtmp = 0:dlB:l_B0-dlB;
    case 'Central'
        xAtmp = dlA/2:dlA:l_A0-dlA/2;
        xBtmp = dlB/2:dlB:l_B0-dlB/2;
    case 'RightBoundary'
        xAtmp = dlA:dlA:l_A0;
        xBtmp = dlB:dlB:l_B0;
end

% l_A and l_B change with each crack increment
l_A = l_A0*ones(q_max+1, 1)+dlB*(0:q_max)';
l_B = l_B0*ones(q_max+1, 1)-dlB*(0:q_max)';

% Matrix where each row spans the entire CLS joint in the (xAB yAB)-Reference
% Frame
xAB = [repmat(xAtmp, q_max+1, 1) repmat(xBtmp, q_max+1, 1)+l_A0];

% Isolate xA and exclude cracked B region elements by setting element index
% to 0; (xA, yA)-Reference Frame
xA = tril(xAB(:, 1:q+q_max), q-1);
xA(tril(ones(size(xA)), q-1)==0) = val;

% Isolate xB and include cracked A region elements; (xB, yB)-Reference Frame
xB = triu(xAB(:, q+1:end)-l_A);
xB(triu(ones(size(xB)))==0) = val;

% xB vector must be adjusted to the x-axis system used in the adhesive
% stress analysis: -l_B <= xB <= 0
xBB = triu(xB-l_B);
xBB(triu(ones(size(xBB)))==0) = NaN;

end