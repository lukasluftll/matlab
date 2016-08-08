function [ref, h, m] = refmap(ls, xgv, ygv, zgv)
% REFMAP Compute reflectivity map from laser scan in grid volume.
%   REF = REFMAP(LS, XGV, YGV, ZGV) uses the laser scan LS to compute the 
%   reflectivity of each voxel in the grid volume defined by the grid 
%   vectors XGV, YGV, ZGV.
%
%   LS is a laserscan object. The sensor pose of the scan is assumed to be 
%   specified with respect to the reflectivity map frame.
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
%   [REF, H, M] = REFMAP(LS, XGV, YGV, ZGV) also returns the voxelmap 
%   objects H and M. 
%   H contains the number of ray remissions for each voxel. 
%   M contains for each voxel the number of rays that traversed the voxel
%   without being reflected.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      ls = laserscan(pcd.azimuth, pcd.elevation, pcd.radius, [1, 100]);
%      ref = refmap(ls, -100:5:100, -100:5:100, -20:5:20)
%
%   See also LASERSCAN, VOXELMAP, REFRAY, DECAYMAP.

% Copyright 2016 Alexander Schaefer
%
% REFMAP implements the counting model proposed by Burgard:
% Wolfram Burgard. Mapping with Known Poses. Lecture Notes on Introduction 
% to Mobile Robotics. University of Freiburg, Germany, 2005.
%http://ais.informatik.uni-freiburg.de/teaching/ss05/robotics/slides/10.pdf

%% Validate input.
% Check number of input arguments.
narginchk(4, 4);

% Check the laserscan.
if ~isa(ls, 'laserscan')
    error('LS must be a laserscan object.')
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

% If the sensor measurement range starts at a positive value, issue a
% warning.
if ls.rlim(1) > 0
    warning(['LS.RLIM(1) > 0, but all no-return ray lengths are ', ...
      'assumed to surpass LS.RLIM(2), not to fall into [0; LS.RLIM(1)].'])
end

%% Preprocess input data.
% Compute the Cartesian ray direction vectors.
ray = cart(ls);

% Compute the indices of the returned rays.
iret = ret(ls);

% Set the length of no-return rays to maximum sensor range.
ray(~iret,:) = ray(~iret,:) * ls.rlim(2);

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
    for i = labindex : numlabs : ls.count
        % Compute the indices of the voxels through which the ray travels.
        [vi, t] = trav(ls.position, ray(i,:), xgv, ygv, zgv);

        % Convert the subscript indices to linear indices.
        vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));

        % For each voxel, sum up the number of hits and misses.
        hw(vi(end)) = hw(vi(end)) + (iret(i) && t(end)==1);
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
