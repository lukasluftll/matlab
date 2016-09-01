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
%   The prior of LF is the Gaussian of the mean distance to the next
%   obstacle inside the map.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);
%      lf = lfmap(pc, 1, -100:5:100, -100:5:100, -20:5:20)
%
%   See also POINTCLOUD, VOXELMAP, LFPC, LFRAY, REFMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(5, 5);

% Check the pointCloud object.
if ~isa(pc, 'pointCloud')
    error('PC must be a pointCloud object.')
end

% Make sure the variance is positive.
validateattributes(sigma, {'numeric'}, {'numel',1, '>',0}, '', 'SIGMA');

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

%% Calculate likelihood field.
% Compute an IxJxK matrix whose rows contain the coordinates of the voxel 
% centers.
[cx, cy, cz] = ndgrid(xgv(1:end-1) + diff(xgv)/2, ...
    ygv(1:end-1) + diff(ygv)/2, zgv(1:end-1) + diff(zgv)/2);
c = [cx(:), cy(:), cz(:)];

% Arrange the points in a 3-column matrix.
points = reshape(pc.Location, pc.Count, 3);

% For each voxel center, compute the distance to the nearest point of the 
% point cloud using multiple workers.
spmd
    % Compute the number of voxels per parallel worker.
    cpw = ceil(size(c,1) / numlabs);

    % Compute the rows of c that this worker considers.
    rc = 1 + (labindex-1)*cpw : min(labindex*cpw, size(c,1));
    
    % Compute the distance to the nearest neighbor.
    [~, dist] = knnsearch(points, c(rc,:));
end

% Merge the results of the individual workers.
dist = reshape(cell2mat(dist(:)), [numel(xgv),numel(ygv),numel(zgv)]-1);

% Compute the likelihood prior.
prior = normpdf(mean(dist(:)), 0, sigma);

% Compute the likelihood field.
lf = voxelmap(normpdf(dist, 0, sigma), xgv, ygv, zgv, prior);

end
