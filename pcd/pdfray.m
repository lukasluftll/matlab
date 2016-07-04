function p = pdfray(origin, ray, lambda, xgv, ygv, zgv)
% PDFRAY Compute probability of Lidar scan given ray decay grid map.
%   P = PDFRAY(ORIGIN, RAY, LAMBDA, XGV, YGV, ZGV) computes the probability
%   of obtaining the Lidar scan defined by ORIGIN and RAY given the decay 
%   map LAMBDA, XGV, YGV, ZGV.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the origins and the 
%   ray vectors of the M rays of the laser scan. Rays that exceed the grid
%   volume are assumed to be no-returns and are attributed the
%   corresponding probability.
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
%      lambda = repmat(magic(5)/100, [1, 1, 5]);
%      gv = 0 : 5; xgv = gv; ygv = gv; zgv = gv;
%      p = pdfray(origin, ray, lambda, xgv, ygv, zgv)
%
%   See also RAYDECAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(6, 6);

% Check if the arguments have the expected sizes.
if size(origin, 1) ~= size(ray, 1)
    error('ORIGIN and RAY must have the same number of rows.')
end
if size(origin, 2) ~= 3 || size(ray, 2) ~= 3
    error('ORIGIN and RAY must have 3 columns.')
end

% Make sure ORIGIN and RAY contain finite values.
if ~all(isfinite(origin)) || ~all(isfinite(ray))
    error('ORIGIN and RAY must not be NaN or Inf.')
end

% Check whether the grid vectors contain enough elements.
if min([numel(xgv), numel(ygv), numel(zgv)]) < 2
    error('Every grid vector must contain at least 2 elements.')
end

% Check whether the grid vectors are ordered.
if any(diff(xgv(:))<=0) || any(diff(ygv(:))<=0) || any(diff(zgv(:))<=0)
    error('Grid vectors must monotonically increase.')
end

% Check whether lambda has the correct size.
if any(size(lambda) ~= [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1])
    error('Size of LAMBDA does not match grid vectors.')
end

%% Compute ray probabilities.
% Determine the number of rays.
nrays = size(origin, 1);

% Loop over all rays.
p = zeros(nrays, 1);
for r = 1 : nrays
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(origin, ray, xgv, ygv, zgv);

    % Compute the lengths of the rays apportioned to each voxel.
    l = diff(t) * norm(ray);

    % Recursively compute the fraction of the rays that arrives at the 
    % beginning of each traversed voxel.
    N = ones(size(l));
    for i = 1 : length(N)-1
        N(i+1) = N(i) * exp(-lambda(vi(i,1),vi(i,2),vi(i,3)) * l(i));
    end
    
    % Compute the probability of obtaining the given ray length depending
    % on whether or not the ray is reflected.
    if t(end) < 1
        p(r) = N(end);
    else
        p(r) = lambda(vi(end,1),vi(end,2),vi(end,3)) * N(end) ...
            * exp(-lambda(vi(end,1),vi(end,2),vi(end,3)) * l(end));
    end
end

end
