function occ = occmap(pc, xgv, ygv, zgv)
% OCCMAP Compute occupancy grid map from point cloud.
%   OCC = OCCMAP(PC, XGV, YGV, ZGV) voxelizes the space according to the 
%   grid vectors XGV, YGV, ZGV and computes a binary occupancy grid map 
%   from the point cloud PC. 
%
%   PC is the pointCloud object.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   OCC is a voxelmap object that contains the binary occupuancy grid. The 
%   value of every voxel that is occupied by at least one point of PC is 
%   true; all other voxels are false.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
%      xgv = min(pcd.x(:)) : 1 : max(pcd.x(:));
%      ygv = min(pcd.y(:)) : 1 : max(pcd.y(:));
%      zgv = min(pcd.z(:)) : 1 : max(pcd.z(:));
%      occ = occmap(pc, xgv, ygv, zgv)
%
%   See also POINTCLOUD, VOXELMAP, LFMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(4, 4);

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

%% Calculate occupancy grid map.
% Create the occupancy grid.
occ = false(numel(xgv)-1, numel(ygv)-1, numel(zgv)-1);

% For each point of the point cloud, compute the corresponding index of the
% voxel grid.
ip = grdidx(reshape(pc.Location, pc.Count, 3), xgv, ygv, zgv);

% Remove rows with zero indices.
ip(any(ip==0, 2),:) = [];

% Set occupied voxels to true.
occ(sub2ind(size(occ), ip(:,1), ip(:,2), ip(:,3))) = true;

% Create voxelmap object.
occ = voxelmap(occ, xgv, ygv, zgv);

end
