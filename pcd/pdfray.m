function L = pdfray(origin, ray, lambda, xgv, ygv, zgv)
% PDFRAY Compute log-likelihood of Lidar measurement given ray decay map.
%   P = PDFRAY(ORIGIN, RAY, LAMBDA, XGV, YGV, ZGV) computes the
%   log-likelihood of obtaining the Lidar ray measurement defined by ORIGIN 
%   and RAY conditioned on the decay map LAMBDA with grid vectors XGV, YGV,
%   ZGV.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the Cartesian 
%   origins and ray vectors of the M measured rays. If all rays originate
%   from the same point, ORIGIN may also be a 1x3 matrix.
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
%   L is an M-element column vector. The value of the m-th element
%   corresponds to the log-likelihood of obtaining the m-th measurement.
%   L is not equal to the log-likelihood of the observation, but it is
%   shifted by an unknown offset.
%
%   Example:
%      origin = [0, 0, 0];
%      ray = [3, 4, 5];
%      lambda = repmat(magic(5)/100, [1, 1, 5]);
%      gv = 0 : 5; xgv = gv; ygv = gv; zgv = gv;
%      p = pdfray(origin, ray, lambda, xgv, ygv, zgv)
%
%   See also NANRAY, RAYDECAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(6, 6);

% Check if the arguments have the expected sizes.
if size(origin, 2) ~= 3 || size(ray, 2) ~= 3
    error('ORIGIN and RAY must have 3 columns.')
end

% If ORIGIN has only one row, expand it to match the row size of RAY.
if size(origin, 1) == 1
    origin = repmat(origin, size(ray, 1), 1);
end

% Make sure all input arguments contain finite values only.
if ~all(isfinite([origin(:); ray(:); lambda(:); xgv(:); ygv(:); zgv(:)]))
    error('All input arguments must not be NaN or Inf.')
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

%% Compute log-likelihood of measurements.
% Determine the number of rays.
nray = size(origin, 1);

% Preallocate the return matrix.
L = zeros(nray, 1);

% Loop over all rays.
parfor i = 1 : nray
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(origin(i,:), ray(i,:), xgv, ygv, zgv);
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(size(lambda), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the length of the ray apportioned to each voxel.
    d = diff(t) * norm(ray(i,:));

    % Compute the log-likelihood of the measurement.
    L(i) = log(lambda(vi(end))) - sum(lambda(vi) .* d);
end

end
