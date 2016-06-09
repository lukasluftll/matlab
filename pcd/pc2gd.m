function [mean, cov] = pc2gd(cloud, res, vol)
% PC2GD Compute Gauss distribution of points in each voxel of grid volume.
%   [MEAN, COV] = PC2GD(CLOUD, RES) divides the axis-aligned volume spanned
%   by the point cloud CLOUD into cubic voxels. Then, for each voxel,
%   the function calculates the mean and the covariance of the point 
%   positions inside the voxel.
%
%   [MEAN, COV] = PC2GD(CLOUD, RES, VOL) does not work on the volume 
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
%   MEAN is a AxBxC matrix that contains the mean position of the points
%   in each voxel. A, B, and C are the counts of voxels in x, y, and z
%   direction. The mean of a voxel that contains no points is set to NaN.
%
%   COV is a AxBxCx3x3 matrix. COV(a,b,c,:,:) is the 3x3 covariance matrix
%   of the voxel with indices a, b, c.
%
%   Example:
%      pc = pcread('teapot.ply');
%      [mean, cov] = pc2gd(pc, 0.1)
%
%   See also POINTCLOUD, NAN.

% Copyright 2016 Alexander Schaefer


end

