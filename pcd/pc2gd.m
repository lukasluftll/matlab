function [mu, sigma] = pc2gd(cloud, res, vol)
% PC2GD Compute Gauss distribution of points in each voxel of grid volume.
%   [MU, SIGMA] = PC2GD(CLOUD, RES) divides the axis-aligned volume 
%   spanned by the point cloud CLOUD into cubic voxels. Then, for each 
%   voxel the function calculates the mean and the covariance of the  
%   point positions inside the voxel.
%
%   [MU, SIGMA] = PC2GD(CLOUD, RES, VOL) does not work on the volume 
%   spanned by CLOUD, but on the axis-aligned cuboid with limits VOL.
%
%   CLOUD is a pointCloud object.
%
%   RES is a scalar that defines the edge length of all voxels that build
%   the grid volume. The voxels are axis-aligned. This means that the edges
%   of the voxels closest to the coordinate axes coincide with the axes.
%   A voxel contains all points [x, y, z]  that satisfy the inequality:
%      (vxmin <= x < vxmax) && (vymin <= y < vymax) && (vzmin <= z < vzmax)
%   with vxmin, vxmax, vymin, vymax, etc. being the limits of the voxel.
%
%   VOL is a 6-element row vector [xmin, ymin, zmin, xmax, ymax, zmax]
%   that describes the limits of the axis-aligned grid volume, including  
%   the minima, excluding the maxima. 
%
%   MU is a AxBxCx3 matrix that contains the mean position of the points
%   in each voxel. A, B, and C are the counts of voxels in x, y, and z
%   direction. MU(a,b,c,:) is the mean point position corresponding to the
%   voxel with indices a, b, c. The mean of a voxel that contains no points 
%   is set to NaN.
%
%   SIGMA is a AxBxCx3x3 matrix. SIGMA(a,b,c,:,:) is the 3x3 covariance 
%   matrix of the voxel with indices a, b, c. The covariance of a voxel 
%   that contains no points is set to NaN.
%
%   Example:
%      pc = pcread('teapot.ply');
%      [mu, sigma] = pc2gd(pc, 0.5)
%
%   See also POINTCLOUD, NAN.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(2, 3);

% If the grid volume is not given, set it to the extent of the point cloud.
if nargin < 3
    vol = reshape([cloud.XLimits; cloud.YLimits; cloud.ZLimits], 1, 6);
    
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

%% Compute the mean and covariance values.
% Calculate the voxel count in each dimension.
voxelcount = ceil(vol(4:6)/res) - floor(vol(1:3)/res);

% Construct the return matrices.
mu = NaN([voxelcount, 3]);
sigma = NaN([voxelcount, 3, 3]);

% Loop over all voxels and detect the numbers of points in each voxel.
for x = 1 : voxelcount(1)
    for y = 1 : voxelcount(2)
        for z = 1 : voxelcount(3)
            % Define the limits of the voxel.
            roi = [x-1, x; y-1, y; z-1, z] * res ...
                + repmat(floor(vol(1:3)'/res) * res, 1, 2);
            
            % Make sure points on the joint face of two voxels are only
            % counted once.
            roi(:,2) = roi(:,2) - eps(roi(:,2));
            
            % Get the points inside the voxel.
            voxelcloud = select(cloud, findPointsInROI(cloud, roi));
            mu(x,y,z,:) = mean(voxelcloud.Location);
            sigma(x,y,z,:,:) = cov(voxelcloud.Location);
        end
    end
end

end
