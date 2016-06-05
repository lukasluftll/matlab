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
%   containing the line parameters that encode where the rays enter and 
%   leave the boxes. 
%   SUPPORT(n,:) + T(1,n)*RAY(n,:) is the point where the n-th ray 
%   intersects with the minimum limits of the n-th box. 
%   SUPPORT(n,:) + T(2,n)*RAY(n,:) is the point where the n-th ray 
%   intersects with the maximum limits of the n-th box; this point does not
%   belong to the box.
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

% Compute the parameters corresponding to the points where the rays enter 
% and leave the box.
t = reshape(t', 3, 2, []);
t(repmat(any(isnan(t), 2), 1, 2)) = NaN;
t = sort(t, 2);
t = reshape([max(t(:,1,:)), min(t(:,2,:))], 2, []);

%% Check for intersection.
% Intersection means the ray travels some distance inside the box.
hit = diff(t) > 0;

% If the ray touches an edge or a corner of the box, this can only be 
% counted as a hit if no coordinate of the intersection represents an upper
% limit of the box.
touch = diff(t) == 0;
contactPoint = support + repmat(t(1,:), 3, 1)' .* ray;
boxmax = box(:, 4:6);
icomp = repmat(touch, 3, 1)';
hit(touch) = ~any(reshape(contactPoint(icomp) == boxmax(icomp), 3, []));

% In case there is no intersection, set t to NaN.
t([~hit; ~hit]) = NaN;

end
