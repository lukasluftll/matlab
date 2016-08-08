function L = lfpc(pc, lf)
% LFPC Compute log-likelihood of point cloud given likelihood field.
%   L = LFPC(PC, LF) computes the likelihood of obtaining the Lidar scan
%   PC given a likelihood field LF.
%
%   PC is a point cloud object; LF is a voxelmap object.
%
%   L is the log-likelihood of obtaining PC.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      pc = pointCloud([pcd.x, pcd.y, pcd.z]);
%      lf = lfmap(pc, 1, xgv, ygv, zgv);
%      l = lfpc(pc, lf)
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

% Check the input argument types.
if ~isa(pc, 'pointCloud')
    error('PC must be a pointCloud object.')
end
if ~isa(lf, 'voxelmap')
    error('MAP must be a voxelmap object.')
end

%% Compute likelihood of point cloud.
% For each point of the point cloud, compute the corresponding index of the
% voxel grid.
ip = grdidx(reshape(pc.Location, pc.Count, 3), lf.xgv, lf.ygv, lf.zgv);

% Remove rows with zero indices.
ip(any(ip==0, 2),:) = [];

% Convert the subscript indices to linear indices.
ip = sub2ind(size(lf.data), ip(:,1), ip(:,2), ip(:,3));

% Compute the log-likelihood of the scan.
L = sum(log(lf.data(ip)));

end
