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
%      data = cat(3, magic(10) / 10^2, magic(10)' / 10^2);
%      alphaplot(data);
%
%   See also: CUBOID, PATCH, CAT.

% Copyright 2016 Alexander Schaefer


narginchk(1, inf)

if ndims(x) ~= 3
    error('X must have exactly 3 dimensions.')
end

if min(x(:)) < 0 || max(x(:)) > 1
    error('X must be element of [0; 1].')
end

if nargin < 2
    vol = [0, 0, 0, size(x)];
end

[x, y, z] = meshgrid(0 : size(x, 1)-1, 0 : size(x, 2)-1, 0 : size(x, 3)-1);
length = (vol(4:6)-vol(1:3)) ./ size(x);
vox = [x(:), y(:), z(:)] .* repmat(length, numel(x), 1) ...
    + repmat(vol(1:3), numel(x), 1);
vox = [vox, vox + repmat(length, size(vox, 1), 1)];

faceAlpha = x(:);
cuboid(vox, 'EdgeColor', 'none', varargin{:}, ...
    'FaceVertexAlphaData', faceAlpha, 'FaceAlpha', 'flat');

end

