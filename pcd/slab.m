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
%   [HIT, T] = SLAB(SUPPORT, RAY, BOX) also returns the 2xN matrix T 
%   containing the line parameters that encode where the rays intersect
%   with the planes that confine the boxes.
%   SUPPORT(n,:) + T(1,n)*RAY(n,:) is the point where the n-th ray enters
%   the n-th box.
%   SUPPORT(n,:) + T(2,n)*RAY(n,:) is the leaving point.
%   If the ray does not intersect with the N-th box, T(:,n) is NaN.
%
%   Example for one ray and one box:
%      support = zeros(1, 3);
%      ray = [1, 1, 1];
%      box = [2, 2, 2, 3, 4, 5];
%      [hit, t] = slab(support, ray, box)
%   
%   Example for multiple rays and multiple boxes:
%      support = zeros(2, 3);
%      ray = [1, 0, 1; 0, 1, 0];
%      box = [1, 0, 0, 2, 1, 1; 10, 10, 10, 15, 20, 11];
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
narginchk(3, 3);

% Check if the input arguments have the expected numbers of columns.
if size(support, 2) ~= 3 || size(ray, 2) ~= 3 || size(box, 2) ~= 6
    error('Input arguments have incorrect row lengths.')
end

% Check if the numbers of rows of the input arguments match.
if size(support, 1) ~= size(ray, 1) || size(support, 1) ~= size(box, 1)
    error('Input argument column lengths must be equal.')
end

% Make sure the box limits are sorted.
reshape(sort([box(1:3), box(4:6)], 2), 6, 1);

%% Compute line parameters of intersection.
% Compute the line parameters of the intersections of the ray and the 
% infinite planes that confine the box.
t = (box - [support, support]) ./ [ray, ray];
t = reshape(t', 3, 2, []);

%% Check for intersection.
% Check if the ray passes through the intersection of at least two 
% lower limit planes.
throughEdge = any([t(1,1,:) == t(2,1,:) & t(1,1,:) ~= t(3,2,:); ...
    t(2,1,:) == t(3,1,:) & t(2,1,:) ~= t(1,2,:); ...
    t(1,1,:) == t(3,1,:) & t(1,1,:) ~= t(2,2,:)]);

% Check if the ray lies inside an upper limit plane.
insideUpper = any(isnan(t(:,2,:)));

% Compute the parameters corresponding to the points where the ray enters 
% and leaves the box.
t(repmat(any(isnan(t), 2), 1, 2)) = NaN;
t = sort(t, 2);
t = reshape([max(t(:,1,:)), min(t(:,2,:))], 2, []);

% The ray intersects with the box if it travels some distance through the 
% box or if it passes through an lower limit edge, but it must not lie 
% inside an upper limit plane.
hit = (diff(t) > 0 | throughEdge) & ~insideUpper;
hit = reshape(hit, 1, []);

% In case there is no intersection, set t to NaN.
t([~hit; ~hit]) = NaN;

end
