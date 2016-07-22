function lf = lfmap(pc, sigma, xgv, ygv, zgv)
% LFMAP Compute likelihood field from point cloud.
%   LF = LFMAP(PC, SIGMA, XGV, YGV, ZGV) transforms the point cloud PC into
%   a voxelized likelihood field using a Gauss kernel with variance SIGMA.
%   The rasterization of the field is defined by the grid vectors XGV, YGV,
%   ZGV.
%
%   PC is the pointCloud object.
%
%   SIGMA is a scalar that defines the variance of the Gaussian kernel.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   LF is a voxelmap object that contains the likelihood of each voxel.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
%      xgv = min(pcd.x(:)) : 1 : max(pcd.x(:));
%      ygv = min(pcd.y(:)) : 1 : max(pcd.y(:));
%      zgv = min(pcd.z(:)) : 1 : max(pcd.z(:));
%      lf = lfmap(pc, 1, xgv, ygv, zgv)
%
%   See also POINTCLOUD, VOXELMAP, LFPC, REFMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(5, 5);

% Make sure the variance is positive.
if sigma <= 0
    error('SIGMA must be positive.')
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

%% Calculate likelihood field.
% Compute an IxJxK matrix whose rows contain the coordinates of the voxel 
% centers.
[cx, cy, cz] = ndgrid(xgv(1:end-1) + diff(xgv)/2, ...
    ygv(1:end-1) + diff(ygv)/2, zgv(1:end-1) + diff(zgv)/2);
c = [cx(:), cy(:), cz(:)];

% Make sure the points are not organized.
points = reshape(pc.Location, pc.Count, 3);

% For each voxel center, compute the distance to the nearest point of the 
% point cloud in parallel using multiple workers.spmd 
spmd
    % Compute the number of voxels per parallel worker.
    cpw = ceil(size(c,1) / numlabs);

    % Compute the rows of c that this worker considers.
    rc = 1 + (labindex-1)*cpw : min(labindex*cpw, size(c,1));
    
    % Compute the distance to the nearest neighbor.
    [~, dist] = knnsearch(points, c(rc,:));
end

% Merge the results of the individual workers.
dist = reshape(cell2mat(dist(:)), [numel(xgv),numel(ygv),numel(zgv)] - 1);

% Compute the likelihood field.
lf = voxelmap(normpdf(dist, 0, sigma), xgv, ygv, zgv);

end
