function [lf] = lfmap(pc, sigma, xgv, ygv, zgv)
% LFMAP Compute likelihood field from point cloud.
%   L = LFMAP(PC, SIGMA, XGV, YGV, ZGV) transforms the point cloud PC into
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
%   LF is a IxJxK matrix that contains the likelihood of each voxel, 
%   where I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
%      xgv = min(pcd.x(:)) : 1 : max(pcd.x(:));
%      ygv = min(pcd.y(:)) : 1 : max(pcd.y(:));
%      zgv = min(pcd.z(:)) : 1 : max(pcd.z(:));
%      lf = lfmap(pc, 1, xgv, ygv, zgv)
%
%   See also DECAYMAP, REFMAP.

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

%% Calculate the likelihood field.
% Determine the size of the voxel grid.
gridsize = [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1];

% Compute a IxJxK matrix whose rows contain the coordinates voxel centers.
[cx, cy, cz] = ndgrid(xgv(1:end-1) + diff(xgv)/2, ...
    ygv(1:end-1) + diff(ygv)/2, zgv(1:end-1) + diff(zgv)/2);
c = [cx(:), cy(:), cz(:)];

% Use multiple workers to loop over all voxels and compute their distances 
% to the nearest point.
spmd
    % Create the distance field for this worker.
    dfw = zeros(gridsize);
    
    for i = labindex : numlabs : size(c, 1)
        [~, dist] = findNearestNeighbors(pc, c(i,:), 1);
        dfw(i) = dist;
    end
end

% Merge the distance fields.
df = zeros(gridsize);
for i = 1 : numel(dfw)
    df = df + dfw{i};
end

% Compute the likelihood field.
lf = normpdf(df, 0, sigma);

end
