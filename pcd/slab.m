function [intersects, t] = slab(origin, ray, box)
% SLAB Compute intersection of ray with axis-aligned box in 3D.
%   INTERSECTS = SLAB(ORIGIN, RAY, BOX) returns whether the infinite ray 
%   characterized by ORIGIN and RAY intersects with the 3D axis-aligned box 
%   described by BOX. 
%
%   ORIGIN and RAY are 3-element column vectors. ORIGIN contains the
%   coordinates of the ray's origin, RAY indicates the direction of the
%   ray. BOX is a 3x2 matrix whose columns contain the minimum and
%   maximum coordinates of the box. The rows correspond to the coordinates.
%
%   [INTERSECTS, T] = SLAB(ORIGIN, RAY, BOX) also returns 
%   a 2-element column vector. ORIGIN + T(1)*RAY is the coordinate of the 
%   point where the ray enters the box, ORIGIN + T(2)*RAY is the point
%   where the ray leaves the box. If the ray does not intersect with the 
%   box, T is [-inf; +inf].
%
%   Example:
%   origin = zeros(3, 1);
%   ray = ones(3, 1);
%   box = [2*ones(3, 1), 3*ones(3, 1)];
%   [intersects, t] = slab(origin, ray, box)
%
%   See also INF.

% Copyright 2016 Alexander Schaefer

% Compute the line parameters of the intersections with the 
% pairwise parallel infinite planes that confine the box. 
t = (box - [origin, origin]) ./ [ray, ray];

% Compute the parameters corresponding to the points where ray enters the
% box and where it leaves it.
t = [max(t(:,1)), min(t(:,2))];

% Check whether some part of the ray remains after computing the entry 
% and leaving point.
intersects = diff(t) >= 0;

end
