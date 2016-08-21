function alphaplot(data, color, xgv, ygv, zgv, varargin)
% ALPHAPLOT Visualize 3D array of scalars using semi-transparent voxels.
%   ALPHAPLOT(DATA) takes the 3-dimensional IxJxK matrix of scalar values 
%   DATA and creates a voxel plot. The element DATA(i,j,k) is represented 
%   by the transparency of the cubic axis-aligned voxel with minimum limits 
%   [i-1,j-1,k-1] and maximum limits [i,j,k].
%   
%   The values of DATA must stay within [0; 1]. DATA(i,j,k) == 0 
%   corresponds to a transparent voxel; a value of 1 corresponds to an 
%   opaque voxel. Values smaller than 0.01 are not drawn.
%
%   ALPHAPLOT(DATA, COLOR) defines the color of each voxel.
%   COLOR can be specified in different ways:
%
%   Single value          - All voxels are plotted in the indexed color 
%                           COLOR.
%
%   3-element row vector  - All voxels are plotted in the true color COLOR.
%
%   IxJxK matrix          - Voxel [i,j,k] is plotted in the indexed color
%                           COLOR(i,j,k).
%
%   3xIxJxK matrix        - Voxel [i,j,k] is plotted in the true color
%                           COLOR(:,i,j,k).
%
%   ALPHAPLOT(DATA, COLOR, XGV, YGV, ZGV) shows no cubic, but cuboid 
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
%      alphaplot(rand(10, 10, 10))
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

% If no color information is given, set it to black.
if nargin < 2
    color = zeros(1, 3);
end

% Check the size of the color matrix.
if ~(numel(color) == 1 ...
        || (ndims(color) == 2 && all(size(color) ~= [1,3])) ...
        || (ndims(color) == 3 && all(size(color) == size(data))) ...
        || (ndims(color) == 4 && all(size(color) == [3, size(data)])))
    error('The size of COLOR is invalid.')
end

% If no volume limits are given, create them.
if nargin < 3
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
[xmin, ymin, zmin] = ndgrid(xgv(1:end-1), ygv(1:end-1), zgv(1:end-1));
[xdiff, ydiff, zdiff] = ndgrid(diff(xgv), diff(ygv), diff(zgv));
vox = [xmin(:), ymin(:), zmin(:)];
vox = [vox, vox + 0.99 * [xdiff(:), ydiff(:), zdiff(:)]];

% Remove voxels with very high transparency.
remove = data(:) < 0.01;
vox(remove,:) = [];
data(remove) = [];

% Define the color values of the voxel faces.
switch ndims(color)
    case 3
        color(remove) = [];
        color = kron(color(:), ones(6,1));
    case 4
        color(:,remove) = [];
        color = kron(reshape(color, 3, []).', ones(6,1));
end

% Define the transparency values of the voxel faces.
faceAlpha = kron(data(:), ones(6,1));

%% Plot voxels.
cuboid(vox, 'EdgeColor', 'none', varargin{:}, ...
    'FaceVertexCData', color, 'FaceColor', 'flat', ...
    'FaceAlpha', 'flat', 'AlphaDataMapping', 'none', ...
    'FaceVertexAlphaData', faceAlpha);
end
