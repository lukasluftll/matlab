function lambda = raydecay(cloud, box)
% RAYDECAY Compute decay rate for rays travelling through box.
%   LAMBDA = RAYDECAY(CLOUD, BOX) uses the pointCloud object CLOUD to 
%   compute mean ray decay rates LAMBDA for each of the boxes specified 
%   by BOX.
%
%   CLOUD is a pointCloud object that contains the readings of a sensor 
%   placed in the origin [0, 0, 0].
%   BOX is a 6xN matrix whose columns describe the limits of each 
%   axis-aligned box [xmin; ymin; zmin; xmax; ymax; zmax].
%   LAMBDA is a N-element row vector containing the mean decay rates for
%   each box.
%
%   Concept of ray decay rate
%   -------------------------
%   The decay rate of a ray emitted by the sensor is a function of the 
%   material through which the ray travels. This function can be 
%   approximated by dividing the space into axis-aligned boxes 
%   and computing the mean decay rate for each box. 
%   The mean decay rate over a box is the number of ray returns from inside
%   the box divided by the sum of the lengths of all rays travelling 
%   through the box:
%
%                 n_returns
%     lambda = ---------------
%               sum(length_i)
%
%   The higher the sum of the ray lengths, the more accurate the
%   approximation of the decay rate.
%
%   RAYDECAY on connected volumes
%   -----------------------------
%   In case RAYDECAY is used on a connected volume of boxes, neighboring 
%   box faces should be REALMIN apart. Otherwise, the neigboring faces of 
%   the two boxes overlap, and returns from inside joint face are counted 
%   twice.
%
% See also POINTCLOUD.

% Copyright 2016 Alexander Schaefer



end
