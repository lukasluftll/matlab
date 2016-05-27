function [intersects, tmin, tmax] = slab(origin, ray, box)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if ray(1) ~= 0
    tmin = (box(1,1) - origin(1)) / ray(1);
    tmax

end

