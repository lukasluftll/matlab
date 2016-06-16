function lambda = raydecay(azimuth, elevation, radius, res, vol)
% RAYDECAY Compute mean decay rate of Lidar rays in voxel.
%   LAMBDA = RAYDECAY(AZIMUTH, ELEVATION, RADIUS, RES) uses the rays 
%   represented in spherical coordinates AZIMUTH, ELEVATION, RADIUS to
%   compute the mean ray decay rate LAMBDA for each of the voxels in the 
%   axis-aligned volume spanned by the rays. The voxels positions are
%   chosen such that no voxel will intersect with the x-y, x-z, or y-z
%   plane. Their edge lengths are given by the scalar RES.
%
%   LAMBDA = RAYDECAY(CLOUD, RES, VOL) computes the mean decay rates of the
%   voxels in a grid volume with resolution RES and extent VOL.
%
%   It is assumed that all rays originate in the origin [0, 0, 0].
%
%   AZIMUTH and ELEVATION are a HEIGHTxWIDTH matrices, where HEIGHT and 
%   WIDTH describe the size of the point cloud. The angle unit is rad.
%
%   RADIUS is a HEIGHTxWIDTH matrix.
%
%   RES is a scalar that defines the edge length of all voxels that build
%   the grid volume. The voxels are axis-aligned. 
%   A voxel contains all points [x, y, z]  that satisfy the inequality:
%      (vxmin <= x < vxmax) && (vymin <= y < vymax) && (vzmin <= z < vzmax)
%   with vxmin, vxmax, vymin, vymax, etc. being the limits of the voxel.
%
%   VOL is a 6-element row vector [xmin, ymin, zmin, xmax, ymax, zmax]
%   that describes the limits of the axis-aligned grid volume, including  
%   the minima, excluding the maxima. 
%
%   LAMBDA is a AxBxC matrix that contains the mean decay rate of each 
%   voxel, where A, C, and C are the counts of voxels in x, y, and z
%   direction. The lambda value of a voxel that has not been visited by any
%   ray is set to NaN.
%
%   Concept of ray decay rate
%   -------------------------
%   The decay rate of a ray emitted by the Lidar sensor is a property of  
%   the material through which the ray travels. This property can be 
%   approximated by dividing the space into voxels and computing the mean 
%   decay rate for each voxel. 
%   The mean decay rate over a voxel is the number of ray returns from 
%   inside the voxel divided by the sum of the lengths of all rays 
%   travelling through the voxel:
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
%      lambda = raydecay(pc.azimuth, pc.elevation, pc.radius, 5)
%
%   See also NAN.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(4, 5);

% Check the spherical coordinate matrices all have the same size.
if any(size(azimuth) ~= size(elevation) | size(azimuth) ~= size(radius))
    error('AZIMUTH, ELEVATION, and RADIUS must all have the same size.')
end

% Check the dimensionality of the spherical coordinate matrices.
if ~(ismatrix(azimuth) && ismatrix(elevation) && ismatrix(radius))
    error('AZIMUTH, ELEVATION, and RADIUS must have exactly 2 dimensions.')
end

% Convert the spherical to Cartesian coordinates.
[x, y, z] = sph2cart(azimuth, elevation, radius);
    
% If the grid volume is not given, set it to the extent of the point cloud.
if nargin < 5
    vol = [min(x(:)), min(y(:)), min(z(:)), ...
        max(x(:)), max(y(:)), max(z(:))];
    
    % Make sure the points lying in the maximum limit planes are included.
    vol(4:6) = vol(4:6) + eps(vol(4:6));
end

% Check the size of the volume vector.
if numel(vol) ~= 6
    error('VOL must have 6 elements.')
end

% Check the volume limits.
if any(diff(reshape(vol', 3, 2), 1, 2) < 0)
    error('Invalid volume limits.')
end

% Check the resolution.
if res <= 0
    error('Resolution must be positive.')
end

%% Sum up ray lengths.
% Construct the matrix that stores the cumulated ray lengths for each
% voxel.
raylength = zeros(ceil(vol(4:6)/res) - floor(vol(1:3)/res));

% Compute the normalized ray direction vectors.
[dirx, diry, dirz] = sph2cart(azimuth, elevation, ones(size(azimuth)));

% Loop over all rays and add the distance they travel through each voxel to
% the matrix that stores the ray lengths.
for i = 1 : numel(azimuth)   
    % Compute the indices of the voxels through which the ray travels.
    [vi, t] = trav([0, 0, 0], [dirx(i), diry(i), dirz(i)], vol, res);
        
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
