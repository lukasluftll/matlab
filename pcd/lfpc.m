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
% Get the coordinates of the points in the cloud as a Nx3 matrix.
p = reshape(pc.Location, pc.Count, 3);

% For each point of the point cloud, compute the corresponding index of the
% voxel grid.
ip = zeros(0, 3);
for i = 1 : size(p,1)
    % Compute the grid index of the point.
    iptmp = [...
        find(map.xgv <= p(i,1) & p(i,1) < map.xgv(end), 1, 'last'), ...
        find(map.ygv <= p(i,2) & p(i,2) < map.ygv(end), 1, 'last'), ...
        find(map.zgv <= p(i,3) & p(i,3) < map.zgv(end), 1, 'last')];
    
    % If the point resides inside the grid, add it to the index matrix.
    ip = [ip; iptmp(repmat(numel(iptmp)==3, size(iptmp)))]; %#ok<AGROW>
end

% Convert the sbuscript indices to linear indices.
ip = sub2ind(size(map.data), ip(:,1), ip(:,2), ip(:,3));

% Compute the log-likelihood of the scan.
L = sum(log(map.data(ip)));

end
