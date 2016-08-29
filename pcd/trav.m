function [i, t] = trav(origin, ray, xgv, ygv, zgv)
% TRAV Ray tracing in voxel grid.
%   I = TRAV(ORIGIN, RAY, XGV, YGV, ZGV) returns the indices of all voxels
%   that a ray traverses on its way from its starting point ORIGIN to its 
%   end point ORIGIN + RAY. The grid is defined by the grid vectors XGV, 
%   YGV, and ZGV.
%
%   If the ray only touches a voxel, i.e. the ray and the voxel have only
%   one point in common, no corresponding traversal is reported.
%
%   ORIGIN is a 3-element row vector that contains the coordinates of
%   the starting point of the ray.
%
%   RAY is a 3-element row vector indicating the direction and the length 
%   of the ray.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   I is an Mx3 matrix whose rows contain the x, y, and z indices of the
%   voxels that the ray traverses, with M being the number of all voxels
%   traversed. I is ordered: The first row corresponds to the voxel in 
%   which the ray starts or where it enters the grid, the last row 
%   corresponds to the voxel where the ray ends or leaves the grid.
%
%   [I, T] = TRAV(ORIGIN, RAY, XGV, YGV, ZGV) also returns the M+1-element 
%   column vector T. It contains the line parameters that encode the 
%   intersections of the ray with the planes that separate the voxels.
%   The first element corresponds to the entry point into the first voxel;
%   in case the ray starts inside a voxel, it is 0. The transition point
%   from the m-th to the m+1st voxel is computed ORIGIN + T(m+1)*RAY.
%   ORIGIN + T(end)*RAY is the point where the ray ends or where it leaves 
%   the grid.
%
%   Example:
%      origin = [5, 1.5, 1];
%      ray = [-10, 0, 0];
%      gv = -2 : 2;
%      [i, t] = trav(origin, ray, gv, gv, gv)
%
%   See also SLAB, NAN.

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
narginchk(5, 5);

% Check if the input arguments are finite.
if ~all(isfinite([origin(:); ray(:); xgv(:); ygv(:); zgv(:)]))
    error('Input arguments must not be NaN or Inf.')
end

% Check if the arguments have the expected sizes.
if numel(origin) ~= 3 || numel(ray) ~= 3
    error('ORIGIN and RAY must have 3 elements.')
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

%% Initialization phase: calculate index of entry point.
% Compute the intersections of the ray with the grid volume.
vol = [xgv(1), ygv(1), zgv(1), xgv(end), ygv(end), zgv(end)];
[hit, tvol] = slab(origin, ray, vol);

% If the ray does not traverse the volume, return immediately.
i = []; 
t = [];
if ~hit || tvol(2) <= 0 || tvol(1) > 1 || diff(tvol) == 0
    return
end

% Compute the line parameters corresponding to the entry point of the ray
% into the grid and to the point where the ray leaves the grid.
tlim = [ max([0, tvol(1)]), min([1, tvol(2)]) ];

% Calculate the index of the voxel corresponding to the starting point.
entry = origin + tlim(1) * ray;
i = [find(letol(xgv(1:end-1), entry(1)), 1, 'last'), ...
    find(letol(ygv(1:end-1), entry(2)), 1, 'last'), ...
    find(letol(zgv(1:end-1), entry(3)), 1, 'last')];

%% Incremental phase: calculate indices of traversed voxels.
% Compute the line parameters of the ray corresponding to the intersections
% with the grid planes.
tx = (xgv(:)-origin(1)) / ray(1);
ty = (ygv(:)-origin(2)) / ray(2);
tz = (zgv(:)-origin(3)) / ray(3);

% Sort the line parameters corresponding to the start point, the end point, 
% and the intersections.
t = unique(sort([0; 1; tx; ty; tz]));

% Remove all intersections that lie outside the grid.
t(t < tlim(1) | t > tlim(2) | ~isfinite(t)) = [];

% For each line parameter compute the index step.
iStep = [ismember(t, tx), ismember(t, ty), ismember(t, tz)] ...
    .* repmat(sign(ray), size(t));

% Sum up the index step to get the voxel indices.
i = cumsum([i; iStep(2:end-1,:)], 1);

end
