function [hit, t] = slab(support, ray, box)
% SLAB Compute intersection of ray with axis-aligned box in 3D.
%   HIT = SLAB(SUPPORT, RAY, BOX) returns whether the infinite ray 
%   characterized by SUPPORT and RAY has at least one point in common with
%   the 3D axis-aligned box described by BOX.
%
%   In case you want to compute the intersection of one ray with one box,
%   SUPPORT and RAY are 3-element column vectors. 
%   SUPPORT contains the coordinates of the ray's point of support.
%   RAY indicates the direction of the ray. 
%   BOX is a 6-element column vector [xmin; ymin; zmin; xmax; ymax; zmax]
%   describing the limits of the box, including the endpoints.
%
%   In case you want to compute N intersections of different rays with 
%   different boxes, SUPPORT and RAY are 3xN matrices.
%   The columns of SUPPORT contain the coordinates of each ray's point of
%   support.
%   The columns of RAY indicate the direction of each ray.
%   BOX is a 6xN matrix whose columns describe the limits of each box.
%   HIT is an N-element logical row vector. 
%
%   [HIT, T] = SLAB(SUPPORT, RAY, BOX) also returns the 2xN matrix T. 
%   SUPPORT(:,N) + T(1,N)*RAY(:,N) is the coordinate of the point where the 
%   N-th ray enters the N-th box;
%   SUPPORT(:,N) + T(2,N)*RAY(:,N) is the point where the ray leaves the 
%   box. 
%   If the N-th ray does not intersect with the N-th box, T(:,N) is NaN.
%
%   SLAB on connected volumes
%   -------------------------
%   In case SLAB is used on a connected volume of boxes, neighboring box
%   faces should be EPS(VERTEX) apart. Otherwise, the neigboring faces of 
%   the two boxes overlap, and rays that travel on the joint face intersect 
%   with both boxes at the same time.
%
%   Example:
%      support = zeros(3, 1);
%      ray = [1; 1; 1];
%      box = [[2; 2; 2]; [3; 3; 3] - eps([3; 3; 3])];
%      [hit, t] = slab(support, ray, box)
%
%   See also NAN, EPS.

% Copyright 2016 Alexander Schaefer
%
% SLAB implements the raycasting algorithm proposed by Smits:
% Brian Smits. Efficiency issues for ray tracing. 
% Journal of Graphics Tools, 3(2):1-14, 1998.

% Compute the line parameters of the intersections of the ray and the 
% infinite planes that confine the box.
t = (box - [support; support]) ./ [ray; ray];
t = sort(reshape(t, 3, 2, []), 2);

% Compute the parameters corresponding to the points where the ray enters 
% and leaves the box.
t = reshape([max(t(:,1,:)), min(t(:,2,:))], 2, []);

% Check whether some part of the ray remains after computing the entry 
% and leaving point.
hit = diff(t) >= 0;

% In case the ray and the box do not intersect, set t to NaN.
t([~hit; ~hit]) = NaN;

end
