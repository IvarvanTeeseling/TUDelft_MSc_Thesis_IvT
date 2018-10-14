function [serr] = DAF_SERR_Effect(x_serr, GI, GII, x_daf, dGI, dGII)

%% Input check

if size(x_serr, 2) ~= 1 || size(x_daf, 2) ~= 1
    error('x_serr and x_daf must be a [Nx1] vector!')
end

if length(x_serr) ~= length(GI) || length(x_serr) ~= length(GII)
    error('x_serr, GI and GII must be the same length!')
end

if length(x_daf) ~= length(dGI) || length(x_daf) ~= length(dGII)
    error('x_daf, dGI and dGII must be the same length!')
end

%% Code

% Use interpolation to find the DAF % effect on the x_serr coordinates
% Note: coordinates outside the x_daf will get 'NaN'
dGI_int     = interp1(x_daf, dGI, x_serr, 'linear');
dGII_int    = interp1(x_daf, dGII, x_serr, 'linear');

% Set all 'NaN' to zero
dGI_int(isnan(dGI_int))     = 0;
dGII_int(isnan(dGII_int))   = 0;

% Include the DAF effect
serr.GI     = GI+GI.*dGI_int;
serr.GII    = GII+GII.*dGII_int;
serr.G      = serr.GI+serr.GII;
serr.MR     = serr.GII./serr.G;

end