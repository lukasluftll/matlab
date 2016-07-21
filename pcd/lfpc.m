function L = lfpc(pc, map)
% LFPC Compute log-likelihood of point cloud given likelihood field.
%   L = LFPC(pc, map) computes the likelihood of obtaining the Lidar scan
%   PC given a likelihood field MAP.
%
%   PC is a point cloud object; MAP is a voxelmap object.
%
%   L is the log-likelihood of obtaining PC.
%
%   Example:
%      pcd = pcdread('data/castle.pcd');
%      pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
%      xgv = min(pcd.x(:)) : max(pcd.x(:));
%      ygv = min(pcd.y(:)) : max(pcd.y(:));
%      zgv = min(pcd.z(:)) : max(pcd.z(:));
%      map = lfmap(pc, 1, xgv, ygv, zgv);
%      l = lfpc(pc, map)
%
%   See also POINTCLOUD, VOXELMAP, LFMAP.

% Copyright 2016 Alexander Schaefer
%
% LFPC implements the likelihood field model as described in:
% S. Thrun, W. Burgard, D. Fox. Probabilistic Robotics. 
% The MIT Press, 2005.

%% Validate input.
% Check number of input arguments.
narginchk(2, 2)

%% Compute likelihood of point cloud.
% For each point of the point cloud, compute the corresponding index of the
% voxel grid.
ip = grdidx(reshape(pc.Location, pc.Count, 3), map.xgv, map.ygv, map.zgv);

% Remove rows with zero indices.
ip(any(ip==0, 2),:) = [];

% Convert the subscript indices to linear indices.
ip = sub2ind(size(map.data), ip(:,1), ip(:,2), ip(:,3));

% Compute the log-likelihood of the scan.
L = sum(log(map.data(ip)));

end
