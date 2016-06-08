function lambda = raydecay(cloud, vol, res)
% RAYDECAY Compute decay rate of different Lidar sensor rays.
%   LAMBDA = RAYDECAY(CLOUD, VOL, RES) uses the pointCloud object CLOUD to 
%   compute the mean ray decay rate LAMBDA for each of the voxels of the 
%   grid volume defined by VOL and RES.
%
%   CLOUD is an organized pointCloud object that contains the readings of 
%   a Lidar sensor. All rays must have originated in the origin [0, 0, 0].
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
%      lambda = raydecay(pc, [pc.XLimits, pc.YLimits, pc.ZLimits)
%
%   See also POINTCLOUD.

% Copyright 2016 Alexander Schaefer

% Compute the intersection of each beam with each box.


end
