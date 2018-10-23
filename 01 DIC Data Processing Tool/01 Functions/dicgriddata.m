function  [intXYZ, extXYZ] = dicgriddata(x, y, z, nx, ny, rec, intmet, gapmet, xc, yc)
% By default, griddata() will interpolate across all coordinates on the
% (x, y) grid within the data convex hull. Any gaps in the data within the
% convex hull will thus automatically be filled and must be restored if
% required.

%% Function

% Mask empty data points (marked by Z = 0) to avoid that they
% locally distort the interpolation
mask = (z==0);
xadj = x(~mask);
yadj = y(~mask);
zadj = z(~mask);

% Create (X, Y) grid based on the overlap region dimensions
X = linspace(rec(1), rec(1)+rec(3), size(x, 2)*nx);
Y = linspace(rec(2), rec(2)+rec(4), size(y, 1)*ny);
[X, Y] = meshgrid(X, Y);

% Pre-fill exterior (x, y) mesh grid
% Note: exterior grid locations with DIC data are marked by extZ = NaN
extZ = NaN(size(X));

% Interpolate the raw data over the (x, y) meshgrid
Z = griddata(xadj, yadj, zadj, X, Y, intmet);

% Find X-coordinates without DIC data (columns filled by NaN)
nancols = any(~isnan(Z), 1);

% Only keep the interior intX-coordinates with DIC data
intX = X(:, nancols);
intY = Y(:, nancols);
intZ = Z(:, nancols);

% Mark exterior extX-coordinates without DIC data with zero 
extZ(:, ~nancols) = 0;

% Find Y-coordinates without DIC data (rows filled by NaN)
nanrows = any(~isnan(Z), 2);

% Only keep the interior intY-coordinates with DIC data
intX(~nanrows, :) = [];
intY(~nanrows, :) = [];
intZ(~nanrows, :) = [];

% Mark exterior extY-coordinates without DIC data with zero 
extZ(~nanrows, :) = 0;

% Empty interpolated Z-data on the location of the bolt if requested
if strcmp(gapmet, 'nofill') && ~isempty(xc) && ~isempty(yc) && length(xc) == length(yc)
    % Find data points inside the boundary
    intXYZcirc = inpolygon(intX, intY, xc, yc);
    
    % Suppress data points inside the boundary
    intZ(intXYZcirc) = NaN;
end

% Aggregate output
intXYZ(:,:,1) = intX;
intXYZ(:,:,2) = intY;
intXYZ(:,:,3) = intZ;

extXYZ(:,:,1) = X;
extXYZ(:,:,2) = Y;
extXYZ(:,:,3) = extZ;

end