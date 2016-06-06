function [i, t] = trav(origin, ray, vol, res)
% TRAV Ray tracing in voxel grid.
%   I = TRAV(ORIGIN, RAY, VOL, RES) returns the indices of all voxels that 
%   a semi-infinite ray traverses on its way from its starting point to the
%   border of the axis-aligned grid.
%
%   ORIGIN is a 3-element row vector that contains the coordinates of
%   the starting point of the ray.
%   RAY is a 3-element row vector indicating the direction of the ray.
%   VOL is a 6-element row vector [xmin, ymin, zmin, xmax, ymax, zmax]
%   that describes the limits of the grid volume, including the minima, 
%   excluding the maxima. A voxel contains all points [x, y, z] that 
%   satisfy the inequality:
%      (vxmin <= x < vxmax) && (vymin <= y < vymax) && (vzmin <= z < vzmax)
%   with vxmin, vxmax, vymin, vymax, etc. being the limits of the voxel.
%   The voxels are axis-aligned. This means that the edges of the voxels 
%   closest to the coordinate axes coincide with the axes.
%   RES is a scalar that defines the edge length of all voxels.
%   I is an Mx3 matrix whose rows contain the x, y, and z indices of the
%   voxels that the ray traverses, with M being the number of all voxels
%   traversed. I is ordered: The first row corresponds to the voxel where 
%   the ray starts, the last row corresponds to the voxel where the ray 
%   leaves the grid.
%
%   [I, T] = TRAV(ORIGIN, RAY, VOL, RES) also returns the M+1-element 
%   column vector T. It contains the line parameters that encode the 
%   intersections of the ray with the planes that separate the voxels.
%   The first element corresponds to the entry point into the first voxel;
%   in case the ray starts inside a voxel, it is NaN. The transition point
%   from the m-th to the m+1st voxel is computed ORIGIN + T(m+1)*RAY.
%   The point where the ray leaves the volume is ORIGIN + T(end)*RAY.
%
%   Example:
%      origin = [5, 2, 2];
%      ray = [-1, 0, 0];
%      vol = [-2, -2, -2, 2, 2, 2];
%      [i, t] = trav(origin, ray, vol, 1)
%
%   See also SLAB, NAN, END.

% Copyright 2016 Alexander Schaefer
%
% TRAV is an advancement on the voxel traversal algorithm proposed by 
% Amanatides and Woo: 
% John Amanatides and Andrew Woo. 
% A Fast Voxel Traversal Algorithm for Ray Tracing. 
% Eurographics 1987, pp. 3-10, 1987.
%
% As an advancement on the paper by Amanatides and Woo, this implementation
% does not limit the index increment from voxel to voxel to one dimension
% (for example [0, 1, 0]), but it also allows steps in two or three 
% dimensions (like [-1, 0, 1] or [1, 1, 1]).

%% Validate input.
% Make sure the user specified enough input arguments.
narginchk(4, 4);

% Check if the arguments have the expected sizes.
inputnames = {'origin', 'ray', 'vol', 'res'};
expectedSize = [1, 3; 1, 3; 1, 6; 1, 1];
for argin = 1 : nargin
    if any(size(eval(inputnames{argin})) ~= expectedSize(argin,:))
        error([upper(inputnames{argin}), ' must be a ', ...
            int2str(expectedSize(argin,1)), 'x', ...
            int2str(expectedSize(argin,2)), ' matrix.'])
    end
end

% Check the volume limits.
if any(diff(reshape(vol', 3, 2), 1, 2) < 0)
    error('Invalid volume limits.')
end

%% Initialization phase: calculate index of point of support.
% Initialize return values.
i = []; t = [];

% Compute the intersections of the ray with the grid volume.
[~, tvol] = slab(origin, ray, vol);

% If the ray does not intersect with the volume, return an empty index 
% matrix.
if tvol(2) < 0
    return
end

% Compute the line parameter of the intersection of the ray with the grid.
if tvol(1) >= 0
    t = tvol(1);
end

% Calculate the index of the starting point.
i(1,:) = floor(origin + t*ray - vol(1:3))' / res + ones(1, 3);

% Compute the bounds of the starting voxel.
voxel = [i(end,:) - ones(1, 3), i(end,:)]' * res + [vol(1:3); vol(1:3)];

%% Incremental phase: calculate indices of traversed voxels.
% Compute the index of the next voxel until the ray leaves the grid.
while true
    % Compute the line parameter of the intersection of the ray with the
    % infinite planes that confine the voxel.
    tvox = (voxel - repmat(origin + t(end)*ray, 1, 2)) ./ [ray; ray];
    
    % Compute the line parameter of the intersection point of the ray with
    % the joint face of the current and the next voxel.
    tvox = max(reshape(tvox, 3, 2), [], 2);
    t = [t; min(tvox)]; %#ok<AGROW>
    
    % Determine the index step into the next voxel.
    iStep = (tvox == t(end)) .* sign(ray);
            
    % Compute the bounds of the next voxel.
    voxel = voxel + [iStep; iStep] * res;
    
    % Check if the next voxel still belongs to the grid volume.
    if ~(all(voxel(1:3) <= vol(4:6)) && all(voxel(4:6) > vol(1:3)))
        return
    end
    
    % Add the index of the voxel to the return matrix.
    i(end+1,:) = i(end,:) + iStep'; %#ok<AGROW>
end
