function [lambda,r,l]=decaymap(o,azimuth,elevation,radius,ret,xgv,ygv,zgv)
% DECAYMAP Compute decay rate of Lidar rays in grid volume.
%   LAMBDA = DECAYMAP(O, AZIMUTH, ELEVATION, RADIUS, RET, XGV, YGV, ZGV)
%   uses the rays represented in spherical coordinates AZIMUTH, ELEVATION, 
%   RADIUS to compute the mean ray decay rate LAMBDA for each voxel in the 
%   grid volume defined by the grid vectors XGV, YGV, ZGV.
%
%   O is a 3-element Cartesian origin vector. All rays originate from this 
%   point.
%
%   AZIMUTH and ELEVATION are HEIGHTxWIDTH matrices, where HEIGHT and WIDTH 
%   describe the size of the point cloud. The unit is rad.
%
%   RADIUS is a HEIGHTxWIDTH matrix that contains the length of the
%   respective ray. For no-return rays, this length must equal the maximum
%   sensor range.
%
%   RET is a HEIGHTxWIDTH logical matrix that indicates whether or not the
%   respective ray was reflected within the sensor range. 
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   LAMBDA is a voxelmap object that contains the mean decay rate of each 
%   voxel. The lambda value of a voxel that has not been visited by any ray
%   is NaN.
%
%   [LAMBDA, R, L] = DECAYMAP(O, AZIMUTH, ELEVATION, RADIUS, XGV, YGV, ZGV)
%   also returns the voxelmap objects R and L. 
%   R contains the number of ray remissions for each voxel. 
%   L contains the cumulated length of all rays that traversed the 
%   respective grid cell.
%
%   Concept of ray decay rate
%   -------------------------
%   The decay rate of a ray emitted by a Lidar sensor is a property of the
%   material through which the ray travels. The mean decay rate over a 
%   voxel is the number of ray returns from inside the voxel divided by the 
%   sum of the lengths of all rays travelling through the voxel:
%
%                  n_returns
%      lambda = -----------------
%                sum(ray_length)
%
%   The higher the sum of the ray lengths, the more accurate the
%   approximation of the decay rate.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      radiusFinite = pcd.radius; radiusFinite(isnan(radiusFinite)) = 130;
%      xgv = min(pcd.x(:)) : 5 : max(pcd.x(:));
%      ygv = min(pcd.y(:)) : 5 : max(pcd.y(:));
%      zgv = min(pcd.z(:)) : 5 : max(pcd.z(:));
%      lambda = decaymap([0,0,0], pcd.azimuth, pcd.elevation, ...
%                      radiusFinite, isfinite(pcd.radius), xgv, ygv, zgv)
%
%   See also VOXELMAP, DECAYRAY, DECAYNANRAY, REFMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(8, 8);

% Check the origin vector.
o = o(:)';
if numel(o) ~= 3 || any(~isfinite(o))
    error('O must be a 3-element real vector.')
end

% Check whether the spherical coordinate matrices and the reflection matrix
% all have the same number of dimensions and the same size.
if ~(ismatrix(azimuth) && ismatrix(elevation) && ismatrix(radius) ...
        && ismatrix(ret))
    error('AZIMUTH, ELEVATION, RADIUS, and RET must be 2D matrices.')
end
if any(size(azimuth) ~= size(elevation) | size(azimuth) ~= size(radius) ...
        | size(azimuth) ~= size(ret))
    error('AZIMUTH, ELEVATION, RADIUS, and RET must have the same size.')
end

% Make sure all input arguments are finite.
if ~all(isfinite([azimuth(:); elevation(:); radius(:); ret(:); ...
        xgv(:); ygv(:); zgv(:)]))
    error('Input arguments must not be NaN or Inf.')
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

%% Preprocess input data.
% Compute the ray direction vectors.
[dirx, diry, dirz] = sph2cart(azimuth, elevation, radius);

% Determine the number of rays.
nray = numel(azimuth);

% Determine the size of the voxel grid.
gridsize = [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1];

%% Sum up ray lengths and returns.
% Use multiple workers to compute the ray length per voxel.
spmd
    % Create return matrices for this worker.
    lw = zeros(gridsize);
    rw = zeros(gridsize);
    
    % For all rays of the worker's share compute the ray length per voxel.
    for i = labindex : numlabs : nray   
        % Compute the indices of the voxels through which the ray travels.
        [vi, t] = trav([0,0,0], [dirx(i),diry(i),dirz(i)], xgv, ygv, zgv);

        % Convert the subscript indices to linear indices.
        vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));

        % Sum up the lengths the rays travel in each voxel.
        lw(vi) = lw(vi) + diff(t) * radius(i);

        % In case of reflection, increment the number of returns.
        rw(vi(end)) = rw(vi(end)) + (ret(i) && t(end)==1);
    end
end

%% Compute ray decay rate.
% Merge the results of all workers.
l = zeros(gridsize);
r = zeros(gridsize);
for i = 1 : numel(lw)
    l = l + lw{i};
    r = r + rw{i};
end

% Compute the decay rate.
lambda = voxelmap(r ./ l, xgv, ygv, zgv);
lambda.data(l == 0) = NaN;

end
