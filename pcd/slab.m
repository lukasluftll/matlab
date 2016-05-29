function [hit, t] = slab(support, ray, box)
% SLAB Compute intersection of ray with axis-aligned box in 3D.
%   HIT = SLAB(SUPPORT, RAY, BOX) returns whether the infinite ray 
%   characterized by SUPPORT and RAY intersects with the 3D axis-aligned 
%   box described by BOX.
%
%   SUPPORT and RAY are 3-element column vectors. 
%   SUPPORT contains the coordinates of the ray's point of support.
%   RAY indicates the direction of the ray. 
%   BOX is a 3x2 matrix whose columns contain the minimum and
%   maximum coordinates of the box. The rows correspond to the coordinates.
%
%   [HIT, T] = SLAB(SUPPORT, RAY, BOX) also returns a 2-element row vector. 
%   SUPPORT + T(1)*RAY is the coordinate of the point where the ray enters 
%   the box, SUPPORT + T(2)*RAY is the point where the ray leaves the box. 
%   If the ray does not intersect with the box, T is NaN.
%
%   Definition of intersection
%   --------------------------
%   "The ray intersects the box" means the ray travels some distance 
%   inside or on the surface of the box. 
%
%   Ray on box surface
%   ------------------
%   In case SLAB is used on a connected volume of voxels, the voxel 
%   surfaces should be REALMIN apart to avoid rays being counted twice 
%   when travelling on the joint face of two voxels.
%
%   Ray touching box corners and edges
%   ----------------------------------
%   In case the ray only touches a corner or an edge, SLAB reports a hit.
%   However, T is zero.
%
%   Example:
%   support = zeros(3, 1);
%   ray = ones(3, 1);
%   box = [2*ones(3, 1), 3*ones(3, 1)];
%   [hit, t] = slab(support, ray, box)
%
%   See also NAN, REALMIN.

% Copyright 2016 Alexander Schaefer
%
% SLAB implements the raycasting algorithm proposed by Smits:
% Brian Smits. Efficiency issues for ray tracing. 
% Journal of Graphics Tools, 3(2):1–14, 1998

% Compute the line parameters of the intersections of the ray and the 
% infinite planes that confine the box. 
t = sort((box - [support, support]) ./ [ray, ray], 2);

% Compute the parameters corresponding to the points where ray enters the
% box and where it leaves it.
t = [max(t(:,1)), min(t(:,2))];

% Check whether some part of the ray remains after computing the entry 
% and leaving point.
hit = diff(t) >= 0;

% In case the ray and the box do not intersect, set t to NaN.
t([~hit, ~hit]) = NaN;

end
