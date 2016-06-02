function i = trav(origin, ray, vol, edge)
% TRAV Ray tracing in voxel grid.
%   I = TRAV(ORIGIN, RAY, VOL, EDGE) returns the indices of all voxels that 
%   a ray traverses on its way from its starting point to the border of the 
%   axis-aligned grid.
%
%   ORIGIN is a 3-element column vector that contains the coordinates of
%   the starting point of the ray.
%   RAY is a 3-element column vector indicating the direction of the ray.
%   VOL is a 6-element column vector [xmin; ymin; zmin; xmax; ymax; zmax]
%   that describes the limits of the grid volume, including the minima, 
%   excluding the maxima. The voxels are axis-aligned. This means that the
%   the edges of the voxels closest to the coordinate axes coincide with 
%   the axes.
%   EDGE is a scalar that defines the edge length of all voxels.
%   I is a Mx3 matrix whose rows contain the x, y, and z indices of the
%   voxels the ray traverses.
%
%   Space spanned by one voxel
%   --------------------------
%   A voxel contains all points [x, y, z] that satisfy the inequality:
%      (vxmin <= x < vxmax) && (vymin <= y < vymax) && (vzmin <= z < vzmax)
%   with vxmin, vxmax, vymin, vymax etc. being the limits of the voxel.
%
%   Example
%      i = trav([-3; 1; 3], [0; 0; 1], [-10; -10; -10; 10; 10; 10], 1)
%
%   See also SLAB.

% Copyright 2016 Alexander Schaefer
%
% TRAV is an advancement on the voxel traversal algorithm proposed by 
% Amanatides and Woo: 
% John Amanatides and Andrew Woo. 
% A Fast Voxel Traversal Algorithm for Ray Tracing. 
% Eurographics 1987, pp. 3-10, 1987.
%
% Instead on limiting the index increment from voxel to voxel to one
% dimension (for example [0 1 0]), this implementation also allows
% steps in two or three dimensions (like [-1 0 1]).

%% Validate input.
% Make sure the user specified enough input arguments.
narginchk(4, 4);

% Check if the arguments have the expected sizes.
inputnames = {'origin', 'ray', 'vol', 'edge'};
expectedSize = [3, 1; 3, 1; 6, 1; 1, 1];
for argin = 1 : nargin
    if any(size(eval(inputnames{argin})) ~= expectedSize(argin,:))
        error([upper(inputnames{argin}), ' must be a ', ...
            int2str(expectedSize(argin,1)), 'x', ...
            int2str(expectedSize(argin,2)), ' matrix.'])
    end
end

% Change the volume vector so it includes the endpoints.
vol(4:6) = vol(4:6) - eps(vol(4:6));

%% Initialization phase: calculate index of point of support.
% Compute the intersections of the ray with the grid volume.
[hit, t] = slab(origin, ray, vol);

% If the ray does not intersect with the volume, return an empty index 
% matrix.
if ~hit || t(2) < 0
    i = [];
    return
end

% If the starting point lies outside the volume, move it towards the 
% volume until it touches its surface.
origin = origin + max(0, t(1)) * ray;

% Calculate the index of the starting point.
i(1,:) = floor((origin - vol(1:3))' / edge + ones(1, 3));

% Compute the bounds of the starting voxel.
voxel = [i(end,:) - ones(1, 3), i(end,:)]' * edge + [vol(1:3); vol(1:3)];

%% Incremental phase: calculate indices of traversed voxels.
% Compute the index of the next voxel until the ray leaves the grid.
while true
    % Compute the line parameter of the intersection of the ray with the
    % infinite planes that confine the voxel.
    t = (voxel - [origin; origin]) ./ [ray; ray];
    
    % Determine the index step into the next voxel.
    tmax = max(reshape(t, 3, 2), [], 2);
    iStep = (tmax == min(tmax)) .* sign(ray);
            
    % Compute the bounds of the voxel we are stepping into.
    voxel = voxel + [iStep; iStep] * edge;
    
    % Check if the new voxel still belongs to the grid volume.
    if ~(all(voxel(1:3) <= vol(4:6)) && all(voxel(4:6) > vol(1:3)))
        return
    end
    
    % Add the index of the voxel to the return matrix.
    i(end+1,:) = i(end,:) + iStep'; %#ok<AGROW>
end
