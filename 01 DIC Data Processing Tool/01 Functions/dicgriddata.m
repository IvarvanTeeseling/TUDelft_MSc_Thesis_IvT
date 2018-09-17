function  [x, y, z, zi, zii] = dicgriddata(X, Y, Z, nx, ny, method, xc, yc)
% By default, griddata() will interpolate across all coordinates on the
% (x, y) grid within the data convex hull. Any gaps in the data within the
% convex hull will thus automatically be filled and must be restored if
% required.

%% Function

% Mask empty data points (marked by Z = 0) to avoid that they
% locally distort the interpolation
mask = (Z==0);
Xadj = X(~mask);
Yadj = Y(~mask);
Zadj = Z(~mask);

% Create the (x, y) meshgrid
x = linspace(min(Xadj(:))*0.99, max(Xadj(:)*0.99), size(X, 2)*nx);
y = linspace(min(Yadj(:))*0.99, max(Yadj(:)*0.99), size(Y, 1)*ny);
[x, y] = meshgrid(x, y);

% Interpolate the raw data over the (x, y) meshgrid
z = griddata(Xadj, Yadj, Zadj, x, y, 'linear');

% Empty the given polygon if requested
if strcmp(method, 'nofill') && ~isempty(xc) && ~isempty(yc) && length(xc) == length(yc)
    % Find data points inside the boundary
    in = inpolygon(x, y, xc, yc);
    % Suppress data points inside the boundary
    z(in) = NaN;
end

% Differentiate z w.r.t. y (1-dimension)
zi = diff(z, 1, 1);
zii = diff(z, 2, 1);

% Remove rows and columns full of NaN

end