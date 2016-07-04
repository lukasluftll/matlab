function p = pdfray(origin, ray, lambda, xgv, ygv, zgv)
% PDFRAY Compute probability of Lidar scan given ray decay grid map.
%   P = PDFRAY(ORIGIN, RAY, LAMBDA, XGV, YGV, ZGV) computes the probability
%   of obtaining the Lidar scan defined by ORIGIN and RAY given the decay 
%   map LAMBDA, XGV, YGV, ZGV.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the origins and the 
%   ray vectors of the M rays of the laser scan.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   LAMBDA is a IxJxK matrix that contains the mean decay rate of each map
%   voxel, where I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
%   The lambda value of a voxel that has not been visited by any ray is 
%   NaN.
%
%   P is a M-element column vector. The value of the m-th element
%   corresponds to the probability of obtaining the m-th measurement.
%
%   Example:
%      origin = [0, 0, 0];
%      ray = [3, 4, 5];
%      lambda = repmat(magic(5), [1, 1, 5]);
%      gv = 1 : 5; xgv = gv; ygv = gv; zgv = gv;
%      p = pdfray(origin, ray, lambda, xgv, ygv, zgv)
%
%   See also RAYDECAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.

%%
% Compute the indices of the grid cells that the ray traverses.
[vi, t] = trav(origin, ray, xgv, ygv, zgv);

vi(t > 1,:) = [];
t(t > 1) = [];
t = t / norm(ray);

% Compute the lengths of the rays apportioned to each voxel.
l = diff(t);

N = ones(size(t));
i = 2;
while t(i) < norm(ray)
    N(i) = N(i-1) * exp(-lambda(i-1)*l(i-1));
end

p = lambda(i) * N(i) * exp(-lambda(i) * t(i));

end
