function lambda = raydecay(cloud, res, vol)
% RAYDECAY Compute decay rate of different Lidar sensor rays.
%   LAMBDA = RAYDECAY(CLOUD, RES) uses the pointCloud object CLOUD to 
%   compute the mean ray decay rate LAMBDA for each of the voxels of the 
%   grid volume spanned by CLOUD with the resolution RES.
%
%   LAMBDA = RAYDECAY(CLOUD, RES, VOL) computes the mean decay rates of the
%   voxels in a grid volume with resolution RES and extent VOL.
%
%   CLOUD is an organized pointCloud object that contains the readings of a
%   Lidar sensor. All rays must have originated in the origin [0, 0, 0].
%   VOL is a 6-element row vector [xmin, ymin, zmin, xmax, ymax, zmax]
%   that describes the limits of the axis-aligned grid volume, including  
%   the minima, excluding the maxima. 
%   RES is a scalar that defines the edge length of all voxels that build
%   the grid volume.
%   A voxel contains all points [x, y, z]  that satisfy the inequality:
%      (vxmin <= x < vxmax) && (vymin <= y < vymax) && (vzmin <= z < vzmax)
%   with vxmin, vxmax, vymin, vymax, etc. being the limits of the voxel.
%   The voxels are axis-aligned. This means that the edges of the voxels 
%   closest to the coordinate axes coincide with the axes.
%   LAMBDA is a XxYxZ matrix that contains the mean decay rate of each 
%   voxel, where X, Y, and Z are the numbers of voxels in x, y, and z
%   direction.
%
%   Concept of ray decay rate
%   -------------------------
%   The decay rate of a ray emitted by the Lidar sensor is a property of  
%   the material through which the ray travels. This property can be 
%   approximated by dividing the space into voxels and computing the mean 
%   decay rate for each voxel. 
%   The mean decay rate over a voxel i is the number of ray returns from 
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
%      pc = pcread('teapot.ply');
%      lambda = raydecay(pc, 0.1)
%
%   See also POINTCLOUD.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(2, 3);

% If the grid volume is not given, use the extent of the point cloud.
if nargin < 3
    vol = reshape([cloud.XLimits; cloud.YLimits; cloud.ZLimits], 1, 6);
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
    error('RES must be positive.')
end

%% Sum up ray lengths.
maxIndex = ceil(vol(4:6)/res) - floor(vol(1:3)/res);
raylength = zeros(maxIndex);
ray = pc2sph(cloud.Location);
[ray(:,:,1), ray(:,:,2), ray(:,:,3)] = sph2cart(...
    ray(:,:,1), ray(:,:,2), ones(size(ray(:,:,3))));

for elevation = 1 : size(ray, 1)
    for azimuth = 1 : size(ray, 2)
        tmpRay = permute(ray(elevation,azimuth,:), [1,3,2]);
        [i, t] = trav(zeros(1, 3), tmpRay, vol, res);
        i = sub2ind(size(raylength), i(:,1), i(:,2), i(:,3));
        raylength(i) = raylength(i) + diff(t);
    end
end

%% Count returns per voxel.
ret = zeros(maxIndex);
for x = 1 : size(raylength, 1)
    for y = 1 : size(raylength, 2)
        for z = 1 : size(raylength, 3)
            roi = [x-1, x; y-1, y; z-1, z]*res ...
                + repmat(floor(vol(1:3)'/res)*res, 1, 2);
            roi(:,2) = roi(:,2) - eps(roi(:,2));
            ret(x,y,z) = ret(x,y,z) + numel(findPointsInROI(cloud, roi));
        end
    end
end

end
