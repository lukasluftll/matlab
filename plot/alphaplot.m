function alphaplot(data, color, xgv, ygv, zgv, varargin)
% ALPHAPLOT Visualize 3D array of scalars using semi-transparent voxels.
%   ALPHAPLOT(DATA) takes the 3-dimensional IxJxK matrix of scalar values 
%   DATA and creates a voxel plot. The element DATA(i,j,k) is represented 
%   by the transparency of the cubic axis-aligned voxel with minimum limits 
%   [i-1,j-1,k-1] and maximum limits [i,j,k].
%   
%   The values of DATA must stay within [0; 1]. A value of 0 corresponds to
%   a transparent voxel; a value of 1 corresponds to an opaque voxel. 
%   Values smaller than 0.01 are not drawn.
%
%   ALPHAPLOT(DATA, COLOR) defines the color of each voxel.
%   COLOR can be specified in different ways:
%
%   Single scalar  - All voxels are plotted in the indexed color COLOR.
%
%   IxJxK matrix   - Voxel [i,j,k] is plotted in color COLOR(i,j,k).
%
%   ALPHAPLOT(DATA, COLOR, XGV, YGV, ZGV) does not show cubic, but cuboid 
%   voxels. XGV, YGV, ZGV are vectors that define the rasterization of the 
%   grid. A voxel with index [i,j,k] contains all points [x,y,z] that 
%   satisfy the inequality:
%
%         (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   ALPHAPLOT(DATA, COLOR, XGV, YGV, ZGV, VARARGIN) plots voxels with the 
%   properties indicated by the name-value pair arguments VARARGIN. 
%   For possible name-value pairs, see the documentation of PATCH.
% 
%   Example:
%      alphaplot(rand(10,10,10))
%
%   See also: CUBOID, PATCH.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user specified the correct number of input arguments.
narginchk(1, inf)

% Check the data matrix is at most 3D.
if ndims(data) > 3
    error('DATA must not have more than 3 dimensions.')
end

% Check the values of the data matrix.
if min(data(:)) < 0 || max(data(:)) > 1
    error('DATA must stay in range [0;1].')
end

% Compute the size of the data matrix.
datasize = [size(data,1), size(data,2), size(data,3)];

% If no color information is given, set it.
if nargin < 2
    color = 0;
end

% If COLOR is a scalar, expand it to match the size of DATA.
if numel(color) == 1
    color = repmat(color, size(data));
end

% Check the size of COLOR.
if ndims(color) > 3
    error('COLOR must not have more than 3 dimensions.')
end
colorsize = [size(color,1), size(color,2), size(color,3)];
if ~(all(colorsize == datasize))
    error('COLOR does not match the size of DATA.')
end

% If no volume limits are given, create them.
if nargin < 3
    xgv = 0 : size(data, 1);
    ygv = 0 : size(data, 2);
    zgv = 0 : size(data, 3);
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

% Check whether the grid vectors match the size of data.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
if ~all(datasize == gridsize)
    error('Grid vectors do not match size of DATA.')
end

%% Plot voxels.
% Compute the limits of the voxels. Make sure the voxel faces do not
% overlap.
[xmin, ymin, zmin] = ndgrid(xgv(1:end-1), ygv(1:end-1), zgv(1:end-1));
[xdiff, ydiff, zdiff] = ndgrid(diff(xgv), diff(ygv), diff(zgv));
vox = [xmin(:), ymin(:), zmin(:)];
vox = [vox, vox + 0.99 * [xdiff(:), ydiff(:), zdiff(:)]];

% Remove voxels with very high transparency.
rem = data(:) < 0.01;
vox(rem,:) = [];
data(rem) = [];
color(rem) = [];

% Define the color values of the voxel faces.
color = kron(color(:), ones(6,1));

% Define the transparency values of the voxel faces.
faceAlpha = kron(data(:), ones(6,1));

%% Plot voxels.
cuboid(vox, 'EdgeColor', 'none', varargin{:}, ...
    'FaceVertexCData', color, 'FaceColor', 'flat', ...
    'FaceAlpha', 'flat', 'AlphaDataMapping', 'none', ...
    'FaceVertexAlphaData', faceAlpha);
end
