function [intersects, inpoint, outpoint] = slab(origin, ray, box)
% SLAB Compute intersection of ray and axis-aligned box in 3D.
%   INTERSECTS = SLAB(ORIGIN, RAY, BOX) returns whether the ray described
%   by ORIGIN and RAY intersects with the 3D axis-aligned box described by
%   BOX. 
%
%   ORIGIN and RAY are 3-element vectors. ORIGIN contains the
%   coordinates of the ray's origin, RAY indicates the direction of the
%   ray. BOX is a 3x2 matrix whose columns contain the minimum and
%   maximum coordinates of the box. The rows correspond to the coordinates
%   x, y, and z.
%
%   [INTERSECTS, INPOINT, OUTPOINT] = SLAB(ORIGIN, RAY, BOX) also returns 
%   the point where the ray enters the box and where it leaves the box.
%   If it does not enter or leave the box, the corresponding value is 
%   set to NaN.
%
%   Example:
%   origin = zeros(3, 1);
%   ray = ones(3, 1);
%   box = [2*ones(3, 1); 3*ones(3, 1)];
%   [intersects, inpoint, outpoint] = slab(origin, ray, box)

% Copyright 2016 Alexander Schaefer

if ray(1) ~= 0
    tmin = (box(1,1) - origin(1)) / ray(1);
    tmax

end

