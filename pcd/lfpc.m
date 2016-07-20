function L = lfpc(pc, map)
% LFPC Compute likelihood of point cloud given likelihood field.
%   L = LFPC(pc, map) computes the likelihood of obtaining the Lidar scan
%   PC given a likelihood field MAP.
%
%   PC is a point cloud object; MAP is a voxelmap object.
%
%   L is the likelihood of obtaining PC.
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
% LFPC implements the likelihood field model like described in:
% S.Thrun, W. Burgard, D. Fox. Probabilistic Robotics. The MIT Press, 2005.

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(2, 2)

%% Compute likelihood of point cloud.
% Get the coordinates of the points in the cloud as a Nx3 matrix.
p = reshape(pc.Location, pc.Count, 3);

% Remove all NaN points.
p(logical(prod(~isfinite(p), 2)),:) = [];

% For each point of the point cloud, compute the corresponding index of the
% voxel grid.
for i = 1 : size(p,1)
    iptmp = [find(map.xgv <= p(i,1), 1, 'last'), ...
        find(map.ygv <= p(i,2), 1, 'last'), ...
        find(map.zgv <= p(i,3), 1, 'last')];
    ip(i, true(1,3) & numel(iptmp)==3) = iptmp(true(1,3) & numel(iptmp)==3);
end
ip(ip(:,1)==0,:) = [];
ip = sub2ind(size(map.data), ip(:,1), ip(:,2), ip(:,3));

% Compute the likelihood of the scan.
L = sum(log(map.data(ip)));

end
