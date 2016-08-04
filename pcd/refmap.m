function [ref, h, m] = refmap(mts, azi, ele, r, ret, xgv, ygv, zgv)
% REFMAP Compute reflectivity map from Lidar rays in grid volume.
%   REF = REFMAP(MTS, AZI, ELE, R, RET, XGV, YGV, ZGV) uses the rays 
%   represented in spherical coordinates azimuth AZI, elevation ELE, and 
%   radius R to compute the reflectivity of each voxel in the grid volume 
%   defined by the grid vectors XGV, YGV, ZGV.
%
%   MTS is an affine3d object that defines the pose of the sensor with
%   respect to the reflectivity map frame.
%
%   AZI and ELE are HEIGHTxWIDTH matrices, where HEIGHT and WIDTH describe 
%   the size of the point cloud. The unit is rad.
%
%   R is a HEIGHTxWIDTH matrix that contains the length of the respective 
%   ray. For no-return rays, this length must equal the maximum sensor 
%   range.
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
%   REF is a voxelmap object that contains the reflectivity of each voxel.
%   The reflectivity is a value in [0; 1]. It indicates the fraction of 
%   rays that are reflected by the voxel compared to all rays that reached 
%   the voxel. If the voxel has not been visited by any ray, its 
%   reflectivity is NaN.
%
%   [REF, H, M] = REFMAP(MTS, AZI, ELE, R, RET, XGV, YGV, ZGV) also returns
%   the voxelmap objects H and M. 
%   H contains the number of ray remissions for each voxel. 
%   M contains for each voxel the number of rays that traversed the voxel
%   without being reflected.
%
%   Example:
%      pc = pcdread('castle.pcd');
%      radiusFinite = pc.radius; radiusFinite(isnan(radiusFinite)) = 130;
%      hgv = -100 : 5 : 100;
%      vgv = -20 : 5 : 20;
%      ref = refmap(affine3d(), pc.azimuth, pc.elevation, radiusFinite, ...
%                   isfinite(pc.radius), hgv, hgv, vgv)
%
%   See also VOXELMAP, REFRAY, REFNANRAY, DECAYMAP.

% Copyright 2016 Alexander Schaefer
%
% REFMAP implements the counting model proposed by Burgard:
% Wolfram Burgard. Mapping with Known Poses. Lecture Notes on Introduction 
% to Mobile Robotics. University of Freiburg, Germany, 2005.
%http://ais.informatik.uni-freiburg.de/teaching/ss05/robotics/slides/10.pdf

%% Validate input.
% Check number of input arguments.
narginchk(8, 8);

% Check the sensor pose.
if ~isa(mts, 'affine3d')
    error('ORG must be an affine3d object.')
end

% Check whether the spherical coordinate matrices and the reflection matrix
% all have the same number of dimensions and the same size.
if ~(ismatrix(azi) && ismatrix(ele) && ismatrix(r) && ismatrix(ret))
    error('AZI, ELE, R, and RET must be 2D matrices.')
end
if any(size(azi)~=size(ele) | size(azi)~=size(r) | size(azi)~=size(ret))
    error('AZI, ELE, R, and RET must have the same size.')
end

% Make sure all input arguments are finite.
if ~all(isfinite([azi(:); ele(:); r(:); ret(:); xgv(:); ygv(:); zgv(:)]))
    error('Input arguments must not be NaN or Inf.')
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

%% Preprocess input data.
% Compute the Cartesian ray direction vectors.
[dirx, diry, dirz] = sph2cart(azi, ele, r);

% Rotate the ray direction vectors according to the sensor orientation.
rot = affine3d([mts.T(1:3,:); 0, 0, 0, 1]);
[dirx, diry, dirz] = transformPointsForward(rot, dirx, diry, dirz);

% Determine the number of rays.
nray = numel(azi);

% Determine the size of the voxel grid.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;

%% Count hits and misses.
% Use multiple workers to compute the number of returns and traversals 
% per voxel for each ray.
spmd
    % Create return matrices.
    hw = zeros(gridsize);
    mw = zeros(gridsize);
    
    % Loop over the worker's share of all rays.
    for i = labindex : numlabs : nray
        % Compute the indices of the voxels through which the ray travels.
        trans = mts.T(4,1:3);
        [vi, t] = trav(trans, [dirx(i), diry(i), dirz(i)], xgv, ygv, zgv);

        % Convert the subscript indices to linear indices.
        vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));

        % For each voxel, sum up the number of hits and misses.
        hw(vi(end)) = hw(vi(end)) + (ret(i) && t(end)==1);
        mw(vi) = mw(vi) + t(2:end)<1;
    end
end

%% Compute reflectivity.
% Merge the hits and misses matrices calculated by the workers.
h = hw{1};
m = mw{1};
for i = 2 : numel(hw)
    h = h + hw{i};
    m = m + mw{i};
end

% For each voxel compute the reflectivity value.
ref = voxelmap(h ./ (h + m), xgv, ygv, zgv);

end
