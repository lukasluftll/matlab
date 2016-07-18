function [ref,h,m] = refmap(azimuth, elevation, radius, ret, xgv, ygv, zgv)
% REFMAP Compute reflectivity map from Lidar rays in grid volume.
%   REF = REFMAP(AZIMUTH, ELEVATION, RADIUS, RET, XGV, YGV, ZGV) uses the 
%   rays represented in spherical coordinates AZIMUTH, ELEVATION, RADIUS to
%   compute the reflectivity of each voxel in the grid volume defined by 
%   the grid vectors XGV, YGV, ZGV.
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
%   REF is a IxJxK matrix that contains the reflectivity of each voxel, 
%   where I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1. The
%   reflectivity is a value in [0; 1]. It indicates the proportion of rays 
%   that are reflected by the voxel. If the voxel has not been visited by
%   any ray, its reflectivity is NaN.
%
%   [REF, H, M] = REFMAP(AZIMUTH, ELEVATION, RADIUS, XGV, YGV, ZGV) also 
%   returns the IxJxK matrices H and M. 
%   H contains the number of ray remissions for each voxel. 
%   M contains for each voxel the number of rays that traversed the voxel
%   without being reflected.
%
%   Example:
%      pc = pcdread('castle.pcd');
%      radiusFinite = pc.radius; radiusFinite(isnan(radiusFinite)) = 130;
%      hgv = -100 : 5 : 100;
%      vgv = -20 : 5 : 20;
%      ref = refmap(pc.azimuth, pc.elevation, radiusFinite, ...
%                   isfinite(pc.radius), hgv, hgv, vgv)
%
%   See also REFRAY, REFNANRAY, DECAYMAP.

% Copyright 2016 Alexander Schaefer
%
% REFMAP implements the counting model proposed by Burgard:
% Wolfram Burgard. Mapping with Known Poses. Lecture Notes on Introduction 
% to Mobile Robotics. University of Freiburg, Germany, 2005.
%http://ais.informatik.uni-freiburg.de/teaching/ss05/robotics/slides/10.pdf

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

%% Count hits and misses.
% Use multiple workers to compute the number of returns and traversals 
% per voxel for each ray.
spmd
    % Create return matrices for this worker.
    hi = zeros(gridsize);
    mi = zeros(gridsize);
    
    % Loop over the worker's share of all rays.
    for i = labindex : numlabs : nray
        % Compute the indices of the voxels through which the ray travels.
        [vi, t] = trav([0,0,0], [dirx(i),diry(i),dirz(i)], xgv, ygv, zgv);

        % Convert the subscript indices to linear indices.
        vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));

        % For each voxel, sum up the number of hits and misses.
        hi(vi(end)) = hi(vi(end)) + (ret(i) && t(end)==1);
        mi(vi) = mi(vi) + t(2:end)<1;
    end
end

%% Compute reflectivity.
% Merge the hits and misses matrices calculated by the workers.
h = sum(reshape(cell2mat(hi(:)), [size(hi{1}), numel(hi)]), 4);
m = sum(reshape(cell2mat(mi(:)), [size(mi{1}), numel(mi)]), 4);

% For each voxel compute the reflectivity value.
ref = h ./ (h + m);

end
