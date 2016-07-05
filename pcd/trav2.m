function [i, t] = trav2(origin, ray, xgv, ygv, zgv)
% TRAV Ray tracing in voxel grid.
%   I = TRAV(ORIGIN, RAY, XGV, YGV, ZGV) returns the indices of all voxels
%   that a ray traverses on its way from its starting point ORIGIN to its 
%   endpoint ORIGIN + RAY. The grid is defined by the grid vectors XGV, 
%   YGV, and ZGV.
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
%   which the ray starts, the last row corresponds to the voxel where the 
%   ray ends or where it leaves the grid.
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
%   TRAV runs on the GPU if at least one input matrix is a gpuArray
%   object. In this case, I is a Nx3 matrix and T is a Nx1 matrix, where
%   N = numel(xgv) + numel(ygv) + numel(zgv). The last N-M rows of I and T 
%   are NaN.   
%
%   Example (execution on CPU):
%      origin = [5, 1.5, 1];
%      ray = [-10, 0, 0];
%      gv = -2 : 2; xgv = gv; ygv = gv; zgv = gv;
%      [i, t] = trav2(origin, ray, xgv, ygv, zgv)
%
%   Example (execution on GPU):
%      origin = gpuArray([5, 1.5, 1]);
%      ray = gpuArray([-10, 0, 0]);
%      gv = gpuArray(-2 : 2); xgv = gv; ygv = gv; zgv = gv;
%      [i, t] = trav(origin, ray, xgv, ygv, zgv)
%
%   See also SLAB, GPUARRAY, NAN.

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

% Check if the arguments have the expected sizes.
if numel(origin) ~= 3 || numel(ray) ~= 3
    error('ORIGIN and RAY must have 3 elements.')
end

% Check the ray vector.
if ~all(isfinite(ray))
    error('Ray values must not be NaN or Inf.')
end

% Check whether the grid vectors contain enough elements.
if min([numel(xgv), numel(ygv), numel(zgv)]) < 2
    error('Every grid vector must contain at least 2 elements.')
end

% Check whether the grid vectors are ordered.
if any(diff(xgv(:))<=0) || any(diff(ygv(:))<=0) || any(diff(zgv(:))<=0)
    error('Grid vectors must monotonically increase.')
end

% Find out whether to execute on the CPU or on the GPU.
gpu = isa(origin,'gpuArray') || isa(ray,'gpuArray') ...
    || isa(xgv,'gpuArray') || isa(ygv,'gpuArray') || isa(zgv,'gpuArray');

%% Initialization phase: calculate index of entry point.
% Compute the size of the voxel grid.
gridsize = [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1];

% Compute the intersections of the ray with the grid volume.
vol = [xgv(1), ygv(1), zgv(1), xgv(end), ygv(end), zgv(end)];
[hit, tvol] = slab(origin, ray, vol);

% If the ray does not intersect with the volume, return immediately.
if ~hit || tvol(2) < 0
    i = []; 
    t = [];
    return
end

% Initialize return matrices. 
i = NaN(sum(gridsize)+3, 3);
t = NaN(sum(gridsize)+3, 1);
if gpu
    i = gpuArray(i);
    t = gpuArray(t);
end

% Define the indices of the end of the payload data of the return matrices.
iEnd = 0;
tEnd = 0;

% Compute the line parameter corresponding to the entry point into the 
% grid.
t(tEnd+1) = max([0, tvol(1)]);
tEnd = tEnd + 1;

% Calculate the index of the voxel corresponding to the starting point. 
entry = origin + t(tEnd)*ray;
iNext = [find(xgv(1:end-1) <= entry(1), 1, 'last'), ...
    find(ygv(1:end-1) <= entry(2), 1, 'last'), ...
    find(zgv(1:end-1) <= entry(3), 1, 'last')];

%% Incremental phase: calculate indices of traversed voxels.
tx = (xgv(:) - origin(1)) / ray(1);
ty = (ygv(:) - origin(2)) / ray(2);
tz = (zgv(:) - origin(3)) / ray(3);

tvol(1) = max([0, tvol(1)]);
tv = min([1, tvol(2)]);
tx(tx <= tvol(1) | tx > tv | ~isfinite(tx)) = [];
ty(ty <= tvol(1) | ty > tv | ~isfinite(ty)) = [];
tz(tz <= tvol(1) | tz > tv | ~isfinite(tz)) = [];

tix = [tx, repmat([sign(ray(1)), 0, 0], size(tx))];
tiy = [ty, repmat([0, sign(ray(2)), 0], size(ty))];
tiz = [tz, repmat([0, 0, sign(ray(3))], size(tz))];

ti = [tix; tiy; tiz];
ti = sortrows(ti);
ti = [tvol(1), iNext; ti];

% multiple index step
i = 2;
while i <= size(ti, 1)
    if ti(i,1) == ti(i-1,1)
        ti(i-1,2:end) = ti(i-1,2:end) + ti(i,2:end);
        ti(i,:) = [];
    else
        i = i + 1;
    end
end

if tvol(2) > 1
    ti(end+1,:) = [1, 0, 0, 0];
end

t = ti(:,1);
i = cumsum(ti(1:end-1,2:end));

end
