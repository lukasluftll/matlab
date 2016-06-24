function lambda = raydecay(azimuth, elevation, radius, xgv, ygv, zgv)
% RAYDECAY Compute decay rate of Lidar rays in grid volume.
%   LAMBDA = RAYDECAY(AZIMUTH, ELEVATION, RADIUS, XGV, YGV, ZGV) uses the
%   rays represented in spherical coordinates AZIMUTH, ELEVATION, RADIUS to
%   compute the mean ray decay rate LAMBDA for each voxel in the grid 
%   volume defined by the grid vectors XGV, YGV, ZGV.
%
%   It is assumed that all rays originate in the origin [0, 0, 0].
%
%   AZIMUTH and ELEVATION are a HEIGHTxWIDTH matrices, where HEIGHT and 
%   WIDTH describe the size of the point cloud. The angle unit is rad.
%
%   RADIUS is a HEIGHTxWIDTH matrix that contains the length of each ray.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%         (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   LAMBDA is a IxJxK matrix that contains the mean decay rate of each 
%   voxel, where I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
%   The lambda value of a voxel that has not been visited by any ray is 
%   NaN.
%
%   Concept of ray decay rate
%   -------------------------
%   The decay rate of a ray emitted by a Lidar sensor is a property of the
%   material through which the ray travels. The mean decay rate over a 
%   voxel is the number of ray returns from inside the voxel divided by the 
%   sum of the lengths of all rays travelling through the voxel:
%
%                  n_returns
%      lambda = -----------------
%                sum(ray_length)
%
%   The higher the sum of the ray lengths, the more accurate the
%   approximation of the decay rate.
%
%   Example:
%      pc = pcdread('castle.pcd');
%      xgv = min(pc.x(:)) : 5 : max(pc.x(:));
%      ygv = min(pc.y(:)) : 5 : max(pc.y(:));
%      zgv = min(pc.z(:)) : 5 : max(pc.z(:));
%      lambda=raydecay(pc.azimuth, pc.elevation, pc.radius, xgv, ygv, zgv)
%
%   See also NAN.

% Copyright 2016 Alexander Schaefer

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

% Check whether the grid vectors contain enough elements.
if min([numel(xgv), numel(ygv), numel(zgv)]) < 2
    error('Every grid vector must contain at least 2 elements.')
end

% Check whether the grid vectors are ordered.
if any(diff(xgv(:))<=0) || any(diff(ygv(:))<=0) || any(diff(zgv(:))<=0)
    error('Grid vectors must monotonically increase.')
end

%% Sum up ray lengths.
% Construct the matrix that stores the cumulated ray lengths for each
% voxel.
raylength = zeros(numel(xgv)-1, numel(ygv)-1, numel(zgv)-1);

% Compute the normalized ray direction vectors.
[dirx, diry, dirz] = sph2cart(azimuth, elevation, ones(size(azimuth)));

% Loop over all rays and add the distance they travel through each voxel to
% the matrix that stores the ray lengths.
for i = 1 : numel(azimuth)   
    % Compute the indices of the voxels through which the ray travels.
    [vi, t] = trav([0, 0, 0], [dirx(i), diry(i), dirz(i)], xgv, ygv, zgv);
        
    % Add the length of the ray that is apportioned to a specific voxel
    % to the cumulated ray length of this voxel.
    vi = sub2ind(size(raylength), vi(:,1), vi(:,2), vi(:,3));
    raylength(vi) = raylength(vi) + diff(t);
end

%% Count returns per voxel.
% Construct the matrix that stores the numbers of returns coming from a
% voxel.
nret = zeros(size(raylength));

% Create a pointCloud object for easy point-per-voxel counting.
cloud = pointCloud([x(:), y(:), z(:)]);

% Loop over all voxels and detect the numbers of points in each voxel.
for x = 1 : size(nret, 1)
    for y = 1 : size(nret, 2)
        for z = 1 : size(nret, 3)
            % Define the limits of the voxel.
            roi = [x-1, x; y-1, y; z-1, z] * res ...
                + repmat(floor(vol(1:3)'/res) * res, 1, 2);

            % Make sure points on the joint face of two voxels are only
            % counted once.
            roi(:,2) = roi(:,2) - eps(roi(:,2));
            
            % Determine the number of points inside the voxel.
            nret(x,y,z) = nret(x,y,z) + numel(findPointsInROI(cloud, roi));
        end
    end
end

%% Compute ray decay rate.
lambda = nret ./ raylength;

end
