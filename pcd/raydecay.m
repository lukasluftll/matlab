function [lambda,r,l] = raydecay(azimuth,elevation,radius,ret,xgv,ygv,zgv)
% RAYDECAY Compute decay rate of Lidar rays in grid volume.
%   LAMBDA = RAYDECAY(AZIMUTH, ELEVATION, RADIUS, RET, XGV, YGV, ZGV) uses
%   the rays represented in spherical coordinates AZIMUTH, ELEVATION, 
%   RADIUS to compute the mean ray decay rate LAMBDA for each voxel in the 
%   grid volume defined by the grid vectors XGV, YGV, ZGV.
%
%   It is assumed that all rays originate in the origin [0, 0, 0].
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
%   LAMBDA is a IxJxK matrix that contains the mean decay rate of each 
%   voxel, where I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
%   The lambda value of a voxel that has not been visited by any ray is 
%   NaN.
%
%   [LAMBDA, R, L] = RAYDECAY(AZIMUTH, ELEVATION, RADIUS, XGV, YGV, ZGV)
%   also returns the IxJxK matrices R and L. 
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
%      pc = pcdread('castle.pcd');
%      radiusFinite = pc.radius; radiusFinite(isnan(radiusFinite)) = 130;
%      xgv = min(pc.x(:)) : 5 : max(pc.x(:));
%      ygv = min(pc.y(:)) : 5 : max(pc.y(:));
%      zgv = min(pc.z(:)) : 5 : max(pc.z(:));
%      lambda=raydecay(pc.azimuth, pc.elevation, radiusFinite, ...
%                      isfinite(pc.radius), xgv, ygv, zgv)
%
%   See also RAYREF, TRAV, SLAB, NAN.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(7, 7);

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
% Construct the matrices that store the cumulated ray lengths and the
% number of returns for each voxel.
l = zeros(gridsize);
r = zeros(gridsize);

% For all rays compute the ray length per voxel and the number of returns
% per voxel.
parfor i = 1 : nray
    % Create temporary return matrices for this ray.
    li = zeros(gridsize);
    ri = zeros(gridsize);
    
    % Compute the indices of the voxels through which the ray travels.
    [vi, t] = trav([0, 0, 0], [dirx(i), diry(i), dirz(i)], xgv, ygv, zgv);
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));
    
    % Add the length of the ray that is apportioned to a specific voxel
    % to the cumulated ray length of this voxel.
    li(vi) = diff(t) * radius(i);

    % If the ray is reflected, increment the number of returns of the voxel 
    % where the ray ends.
    ri(vi(end)) = ret(i) && t(end)==1;
    
    % Add the temporary values to the return matrices.
    l = l + li;
    r = r + ri;
end

%% Compute ray decay rate.
lambda = r ./ l;
lambda(l == 0) = NaN;

end
