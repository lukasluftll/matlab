function alphaplot(data, xgv, ygv, zgv, varargin)
% ALPHAPLOT Visualize 3D array of scalars using semi-transparent voxels.
%   ALPHAPLOT(DATA) takes the 3-dimensional array of scalar values DATA and
%   creates a voxel plot. The element DATA(x,y,z) is represented by the 
%   transparency of the cubic axis-aligned voxel with minimum limits 
%   [x-1, y-1, z-1] and maximum limits [x, y, z].
%   
%   The values of DATA must stay within [0; 1]. DATA(x,y,z) == 0 
%   corresponds to a transparent voxel; a value of 1 corresponds to an 
%   opaque voxel. Values smaller than 0.01 are not drawn.
%
%   ALPHAPLOT(DATA, XGV, YGV, ZGV) shows no cubic, but cuboid voxels.
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%         (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   ALPHAPLOT(DATA, XGV, YGV, ZGV, VARARGIN) plots voxels with the 
%   properties indicated by the name-value pair arguments VARARGIN. 
%   For possible name-value pairs, see the documentation of PATCH.
% 
%   Example:
%      data = [0.1, 0.2, 0.3; 0, 0, 0; 0, 0, 0];
%      data = cat(3, data, zeros(3), 0.5 * ones(3));
%      alphaplot(data);
%
%   See also: CUBOID, PATCH, CAT.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user specified the correct number of input arguments.
narginchk(1, inf)

% Check if the data matrix has 3 dimensions.
if ndims(data) ~= 3
    error('DATA must have exactly 3 dimensions.')
end

% Make sure the values of DATA are finite and stay in [0; 1].
if min(data(:)) < 0 || max(data(:)) > 1
    error('X must be element of [0; 1].')
end

% If no volume limits are given, create them.
if nargin < 2
    xgv = 0 : size(data, 1);
    ygv = 0 : size(data, 2);
    zgv = 0 : size(data, 3);
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

% Check whether DATA has the correct size.
if any(size(data) ~= [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1])
    error('Size of DATA does not match grid vectors.')
end

%% Plot voxels.
% Compute the limits of the voxels. Make sure the voxel faces do not
% overlap.
[xmin, ymin, zmin] = meshgrid(xgv(1:end-1), ygv(1:end-1), zgv(1:end-1));
[xdiff, ydiff, zdiff] = meshgrid(diff(xgv), diff(ygv), diff(zgv));
vox = [xmin(:), ymin(:), zmin(:)];
vox = [vox, vox + 0.99 * [xdiff(:), ydiff(:), zdiff(:)]];

% Remove voxels with very high transparency.
remove = data(:) < 0.01;
vox(remove,:) = [];
data(remove) = [];

% Define the transparency values of the voxel faces.
faceAlpha = kron(data(:), ones(6, 1));

%% Plot voxels.
cuboid(vox, 'FaceColor', 'blue', 'EdgeColor', 'none', varargin{:}, ...
    'FaceVertexAlphaData', faceAlpha, 'FaceAlpha', 'flat');

end
