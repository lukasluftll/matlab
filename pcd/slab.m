function [hit, t] = slab(support, ray, box)
% SLAB Compute intersection of ray with axis-aligned box in 3D.
%   HIT = SLAB(SUPPORT, RAY, BOX) returns whether the infinite ray 
%   characterized by SUPPORT and RAY has at least one point in common with
%   the 3D axis-aligned box described by BOX.
%
%   In case you want to compute the intersection of one ray with one box,
%   SUPPORT and RAY are 3-element row vectors. 
%   SUPPORT contains the coordinates of the ray's point of support.
%   RAY indicates the direction of the ray. 
%   BOX is a 6-element row vector [xmin, ymin, zmin, xmax, ymax, zmax]
%   describing the limits of the box, including the minima, excluding the
%   maxima. A point [x, y, z] lies inside the box if it satisfies the 
%   inequalities:
%      (xmin <= x < ymax) && (ymin <= y < ymax) && (zmin <= z < zmax)
%
%   In case you want to compute N intersections of different rays with 
%   different boxes, SUPPORT and RAY are Nx3 matrices.
%   The rows of SUPPORT contain the coordinates of each ray's point of
%   support.
%   The rows of RAY indicate the direction of each ray.
%   BOX is a Nx6 matrix whose rows describe the limits of each box.
%   HIT is an N-element logical row vector. 
%
%   [HIT, T] = SLAB(SUPPORT, RAY, BOX) also returns the Nx2 matrix T. 
%   SUPPORT(N,:) + T(N,1)*RAY(N,:) is the coordinate of the point where the 
%   N-th ray enters the N-th box;
%   SUPPORT(N,:) + T(N,2)*RAY(N,:) is the point where the ray leaves the 
%   box. 
%   If the N-th ray does not intersect with the N-th box, T(N,:) is NaN.
%
%   Example for one ray and one box:
%      support = zeros(1, 3);
%      ray = [1, 1, 1];
%      box = [2, 2, 2, 3, 4, 5];
%      [hit, t] = slab(support, ray, box)
%   
%   Example for multiple rays and multiple boxes:
%      support = zeros(2, 3);
%      ray = [1, 1, 1; 0, 1, 0];
%      box = [2, 2, 2, 3, 4, 5; -10, -10, -10, 15, 20, 11];
%      [hit, t] = slab(support, ray, box)
%
%   See also NAN, TRAV.

% Copyright 2016 Alexander Schaefer
%
% SLAB implements the raycasting algorithm proposed by Smits:
% Brian Smits. Efficiency issues for ray tracing.
% Journal of Graphics Tools, 3(2):1-14, 1998.

%% Validate input.
% Check the number of input arguments.


%% Compute intersection.
% Compute the line parameters of the intersections of the ray and the 
% infinite planes that confine the box.
t = (box - [support, support]) ./ [ray, ray];

% Shorten the ray parameters by an infinite amount to account for the fact
% that the rays have already left the box when they intersect with
% the planes describing the maximum coordinates of the box.
epsilon = -eps(t(:,4:6)) .* sign(ray);
epsilon(isnan(epsilon)) = 0;
t(:,4:6) = t(:,4:6) + epsilon;

% Compute the parameters corresponding to the points where the ray enters 
% and leaves the box.
t = sort(reshape(t', 3, 2, []), 2);
t = reshape([max(t(:,1,:)), min(t(:,2,:))], 2, []);

% Check whether some part of the ray remains after computing the entry 
% and leaving point.
hit = diff(t) >= 0;

% In case the ray and the box do not intersect, set t to NaN.
t([~hit; ~hit]) = NaN;

end
