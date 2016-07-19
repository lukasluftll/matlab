function alphaplot(data, xgv, ygv, zgv, varargin)
% ALPHAPLOT Visualize 3D array of scalars using semi-transparent voxels.
%   ALPHAPLOT(DATA) takes the 3-dimensional array of scalar values DATA and
%   creates a voxel plot. The element DATA(x,y,z) is represented by the 
%   transparency of the cubic axis-aligned voxel with minimum limits 
%   [x-1, y-1, z-1] and maximum limits [x, y, z].
%   
%   The values of DATA must stay within [0; 1]. DATA(x,y,z) == 0 
%   corresponds to a transparent voxel; a value of 1 corresponds to an 
%   opaque voxel.
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
%      data = cat(3, magic(10) * 1e-3, magic(10)' * 1e-3);
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

%% Compute voxel limits.
% Compute an Nx3 matrix whose columns contain the indices of all voxels.
[xmin, ymin, zmin] = ndgrid(...
    1 : size(data, 1), 1 : size(data, 2), 1 : size(data, 3));
i = [xmin(:), ymin(:), zmin(:)];

% Remove the voxels with very high transparency.
remove = data(:) < 0.01;
i(remove,:) = [];
data(remove) = [];

% Compute the limits of the voxels. Make sure the voxel faces do not
% overlap.
minvox = [xgv(i(:,1)); ygv(i(:,2)); zgv(i(:,3))]; ...
maxvox = [xgv(i(:,1)+1); ygv(i(:,2)+1); zgv(i(:,3)+1)];
vox = [minvox; minvox + 0.99*(maxvox-minvox)]';

% Define the transparency values of the voxel faces.
faceAlpha = kron(data(:), ones(6, 1));

%% Plot voxels.
cuboid(vox, 'FaceColor', 'blue', 'EdgeColor', 'none', varargin{:}, ...
    'FaceVertexAlphaData', faceAlpha, 'FaceAlpha', 'flat');

end
