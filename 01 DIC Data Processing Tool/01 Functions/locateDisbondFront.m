function [zb, bA, zp] = locateDisbondFront(x, y, z, search, gap, xc, yc)

% Pre-allocate memory
zb = zeros(size(y));
zp = zeros(size(y));

for i = 1:size(z, 2)
    switch search
        case 'findvalley'
            
            % Get the peak values
            [maxval, maxpos] = findpeaks(z(:,i));
            [minval, minpos] = findpeaks(z(:,i).^-1);
            
            % Merge the found peaks together
            merged = [maxval maxpos ones(size(maxval)) ; ...
                minval.^-1 minpos zeros(size(minval))];
            merged = sortrows(merged, 2);
            
            % Find the location of the maximum difference between a maximum
            % and its previous minimum
            if merged(1,3) == 0
                dz = diff(merged(:, 1));
                ybpos = merged(dz == max(dz), 2);
                zbval = merged(dz == max(dz), 1);
            elseif merged(1,3) == 1
                dz = diff(merged(2:end, 1));
                ybpos = merged(find(dz == max(dz))+1, 2);
                zbval = merged(find(dz == max(dz))+1, 1);
            end
            
%             h =  findobj('type','figure');
%             n = length(h);
%             figure(n+1)
%             hold on
%             plot(y(:,i), z(:,i));
%             plot(y(merged(:,2)), merged(:,1), 'ro')
%             plot(y(ybpos,2), zbval, 'gs')
%             hold off
            
            % Write valley z value to all disbonded y positions
            zb(1:ybpos, i) = 1;
            
            % Save the peak values on the z grid
            zp(merged(:,2)) = merged(:,1);
            
        case 'dz/dy'
        case 'd^2z/dy^2'
    end
end

% Empty the given polygon if requested
if strcmp(gap, 'nofill') && ~isempty(xc) && ~isempty(yc) && length(xc) == length(yc)
    % Find data points inside the boundary
    in = inpolygon(x, y, xc, yc);
    % Suppress data points inside the boundary
    zb(in) = 0;
    zp(in) = 0;
end

% Disbond area
bA = 1;

end