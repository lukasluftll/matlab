function alphaplot(data, vol, varargin)
% ALPHAPLOT Visualize 3D array of scalars using semi-transparent voxels.
%   ALPHAPLOT(DATA) takes the 3-dimensional array of scalar values DATA and
%   creates a voxel plot. The element DATA(x,y,z) is represented by the 
%   transparency of the cubic axis-aligned voxel with minimum limits 
%   [x-1, y-1, z-1] and maximum limits [x, y, z]. 
%   
%   The values of DATA must stay within [0; 1]. DATA(x,y,z) == 0 
%   corresponds to a completely transparent voxel; a value of 1 corresponds
%   to a completely opaque voxel.
%
%   ALPHAPLOT(DATA, VOL) divides the axis-aligned volume VOL into cuboid 
%   voxels and uses this volume to visualize DATA. VOL is a 6-element 
%   vector that contains the minimum and maximum limits of the volume: 
%   [xmin, ymin, zmin, xmax, ymax, zmax].
%
%   ALPHAPLOT(DATA, VOL, VARARGIN) plots voxels with the properties 
%   indicated by the name-value pair arguments VARARGIN. 
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

% If no volume limits are given, use cubic voxels with edge length 1.
if nargin < 2
    vol = [0, 0, 0, size(data)];
end

% Check the volume limits.
if any((vol(4:6) - vol(1:3)) <= 0)
    error('Invalid argument VOL: VOL(1:3) >= VOL(4:6).')
end

%% Compute voxel limits.
% Compute a Nx3 matrix whose columns contain the indices of all voxels.
[xmin, ymin, zmin] = meshgrid(...
    1 : size(data, 1), 1 : size(data, 2), 1 : size(data, 3));
i = [xmin(:), ymin(:), zmin(:)];

% Compute the edge lengths of the voxels.
edgelength = (vol(4:6)-vol(1:3)) ./ size(data);

% Compute the minimum limits of the voxels.
vox = repmat(vol(1:3), numel(xmin), 1) ...
    + (i-1) .* repmat(edgelength, numel(xmin), 1);

% Add the maximum limits.
vox = [vox, vox + repmat(edgelength, size(vox, 1), 1)];

% Change the order of the data elements to the order of the voxels.
faceAlpha = data(sub2ind(size(data), i(:,1), i(:,2), i(:,3)));
faceAlpha = kron(faceAlpha, ones(6, 1));

%% Plot voxels.
cuboid(vox, 'FaceColor', 'blue', 'EdgeColor', 'none', varargin{:}, ...
    'FaceVertexAlphaData', faceAlpha, 'FaceAlpha', 'flat');

end
