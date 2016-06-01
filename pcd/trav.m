function i = trav(origin, ray, vol, edge)
% TRAV Ray tracing in voxel grid.
%   I = TRAV(ORIGIN, RAY, VOL, EDGE) returns the indices of all voxels that 
%   a ray traverses on its way from its starting point to the border of the 
%   axis-aligned grid.
%
%   ORIGIN is a 3-element column vector that contains the coordinates of
%   the ray's starting point.
%   RAY is a 3-element column vector indicating the direction of the ray.
%   VOL is a 6-element column vector [xmin; ymin; zmin; xmax; ymax; zmax]
%   that describes the limits of the grid volume, including the starting 
%   points, excluding the endpoints. The voxels are axis-aligned. This 
%   means that the edges of the voxels closest to the coordinate axes 
%   coincide with the axes.
%   EDGE is a scalar that defines the edge length of the voxels.
%   I is a Mx3 matrix whose rows contain the x, y, and z indices of the
%   voxels the ray traverses.
%
%   Space spanned by one voxel
%   --------------------------
%   A voxel contains all points [x, y, z] that satisfy all of the following
%   inequalities:
%      xmin <= x < xmax 
%      ymin <= y < ymax
%      zmin <= z < zmax
%
%   Example
%      i = trav([-3; 1; 3], [0; 0; 1], [-10; -10; -10; 10; 10; 10], 1)
%
%   See also SLAB.

% Copyright 2016 Alexander Schaefer
%
% TRAV implements the voxel traversal algorithm proposed by Amanatides and 
% Woo:
% John Amanatides and Andrew Woo. 
% A Fast Voxel Traversal Algorithm for Ray Tracing.
% Eurographics 1987, pp. 3-10, 1987.

%% Validate input.
% Make sure the user specified enough input arguments.
narginchk(4, 4);

% Check if the arguments have the expected sizes.
expectedSize = [3, 1; 3, 1; 6, 1; 1, 1];
for argin = 1 : nargin
    if size(eval(inputname(argin))) ~= expectedSize(argin,:)
        error([inputname(argin), ' must be a ', expectedSize(argin,1), ...
            'x', expectedSize(argin,2), ' matrix.'])
    end
end

%% Initialization phase: calculate index of point of support.
% Compute the intersections of the ray with the volume.
[hit, t] = slab(origin, ray, vol);

% If the ray does not intersect with the volume, return an empty index 
% matrix.
if ~hit || t(2) <= 0
    i = [];
    return
end

% If the starting point lies outside the volume, move it towards the 
% volume until it touches its surface.
origin = origin + t([t(1)>0, false]) * ray;

% Calculate the index of the starting point.
i(:,1) = floor((origin - vol(1:3)) ./ edge + ones(3, 1));

%% Incremental phase: calculate indices of traversed voxels.
% Compute the next index until the ray leaves the grid.
while all(voxel(4:6) >= vol(1:3)) || all(voxel(1:3) <= vol(4:6))
    % Compute the bounds of the current voxel.
    voxel = [i(:,end)-ones(3, 1); i(:,end)] .* edge + vol([1:3, 1:3]);
    if voxel > vol
        return
    end
    % Compute the line parameters of the intersections of the ray and the
% infinite planes that confine the voxel.
    t = (voxel - [origin; origin]) ./ [ray; ray];
    
    [~, m] = max(t);
    mplus = [0; 0; 0];
    mplus(m) = 1;
    i = [i, i(:,end) + mplus]; %#ok<AGROW>
    
end
