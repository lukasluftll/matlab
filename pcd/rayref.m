function [ref, h, m] = rayref(azimuth, elevation, radius, xgv, ygv, zgv)
% RAYREF Compute reflectivity map from Lidar rays in grid volume.
%   REF = RAYREF(AZIMUTH, ELEVATION, RADIUS, XGV, YGV, ZGV) uses the rays 
%   represented in spherical coordinates AZIMUTH, ELEVATION, RADIUS to
%   compute the reflectivity of each voxel in the grid volume defined by 
%   the grid vectors XGV, YGV, ZGV.
%
%   It is assumed that all rays originate in the origin [0, 0, 0].
%
%   AZIMUTH and ELEVATION are HEIGHTxWIDTH matrices, where HEIGHT and 
%   WIDTH describe the size of the point cloud. The angle unit is rad.
%
%   RADIUS is a HEIGHTxWIDTH matrix that contains the length of each ray.
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
%   [REF, H, M] = RAYREF(AZIMUTH, ELEVATION, RADIUS, XGV, YGV, ZGV) also 
%   returns the IxJxK matrices H and M. 
%   H contains the number of ray remissions for each voxel. 
%   M contains for each voxel the number of rays that traversed the voxel
%   without being reflected.
%
%   Example:
%      pc = pcdread('castle.pcd');
%      pc.radius(~isfinite(pc.radius)) = 130;
%      xgv = min(pc.x(:)) : 5 : max(pc.x(:));
%      ygv = min(pc.y(:)) : 5 : max(pc.y(:));
%      zgv = min(pc.z(:)) : 5 : max(pc.z(:));
%      ref = rayref(pc.azimuth, pc.elevation, pc.radius, xgv, ygv, zgv)
%
%   See also RAYDECAY, TRAV, SLAB, NAN.

% Copyright 2016 Alexander Schaefer
%
% RAYREF implements the counting model proposed by Burgard:
% Wolfram Burgard. Mapping with Known Poses. Lecture Notes on Introduction 
% to Mobile Robotics. University of Freiburg, Germany, 2005.
%http://ais.informatik.uni-freiburg.de/teaching/ss05/robotics/slides/10.pdf

%% Validate input.
% Check number of input arguments.
narginchk(6, 6);

% Check whether the spherical coordinate matrices all have the same size.
if any(size(azimuth) ~= size(elevation) | size(azimuth) ~= size(radius))
    error('AZIMUTH, ELEVATION, and RADIUS must all have the same size.')
end

% Check the dimensionality of the spherical coordinate matrices.
if ~(ismatrix(azimuth) && ismatrix(elevation) && ismatrix(radius))
    error('AZIMUTH, ELEVATION, and RADIUS must have exactly 2 dimensions.')
end

% Make sure all input arguments are finite.
if ~all(isfinite([azimuth(:);elevation(:);radius(:);xgv(:);ygv(:);zgv(:)]))
    error('Input arguments must not be NaN or Inf.')
end

% Check whether the grid vectors contain enough elements.
gvchk(xgv, ygv, zgv);

%% Preprocess input data.
% Compute the ray direction vectors.
[dirx, diry, dirz] = sph2cart(azimuth, elevation, radius);

% Determine the number of rays.
nray = numel(azimuth);

% Determine the size of the voxel grid.
gridsize = [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1];

%% Sum up ray lengths and returns.
% Construct the matrices that store the cumulated ray lengths and the
% number of returns for each voxel.
m = zeros(gridsize);
h = zeros(gridsize);

% For all rays compute the ray length per voxel and the number of returns
% per voxel.
parfor i = 1 : nray
    % Create temporary return matrices for this ray.
    li = zeros(gridsize);
    ri = zeros(gridsize);
    
    % Compute the indices of the voxels through which the ray travels.
    [vi, t] = trav([0, 0, 0], [dirx(i), diry(i), dirz(i)], xgv, ygv, zgv);
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));
    
    % Add the length of the ray that is apportioned to a specific voxel
    % to the cumulated ray length of this voxel.
    li(vi) = diff(t) * radius(i);

    % Increment the number of returns of the voxel where the ray ends.
    ri(vi(end)) = t(end)==1;
    
    % Add the temporary values to the return matrices.
    m = m + li;
    h = h + ri;
end

%% Compute ray decay rate.
ref = h ./ m;
ref(m == 0) = NaN;

end
