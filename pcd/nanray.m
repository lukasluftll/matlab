function p = nanray(origin, ray, rlim, lambda, xgv, ygv, zgv)
% NANRAY Compute probability of NaN Lidar measurement from decay map.
%   P = NANRAY(ORIGIN, RAY, RLIM, LAMBDA, XGV, YGV, ZGV) computes the
%   probability of obtaining the measurement NaN from a Lidar sensor that 
%   sends a ray from ORIGIN in direction RAY through a ray decay voxel 
%   grid defined by LAMBDA with grid vectors XGV, YGV, ZGV.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the Cartesian 
%   origins and ray vectors of the M measured rays. If all rays originate
%   from the same point, ORIGIN may also be a 1x3 matrix. The length of RAY 
%   is not considered.
%
%   RLIM is a 2-element vector that defines the minimum and the maximum
%   radius detected by the Lidar sensor. All other radii are assumed to
%   result in an NaN measurement.
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
%   P is an M-element column vector. The value of the m-th element
%   corresponds to the log-likelihood of obtaining NaN for the m-th 
%   measurement.
%
%   Example:
%      origin = [0, 0, 0];
%      ray = [3, 4, 5];
%      rlim = [1; 100];
%      lambda = repmat(magic(5)/100, [1, 1, 5]);
%      gv = 0 : 5; xgv = gv; ygv = gv; zgv = gv;
%      p = nanray(origin, ray, rlim, lambda, xgv, ygv, zgv)
%
%   See also PDFRAY, RAYDECAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(7, 7);

% Check if the arguments have the expected sizes.
if size(origin, 2) ~= 3 || size(ray, 2) ~= 3
    error('ORIGIN and RAY must have 3 columns.')
end

% If ORIGIN has only one row, expand it to match the row size of RAY.
if size(origin, 1) == 1
    origin = repmat(origin, size(ray, 1), 1);
end

% Make sure all input arguments contain finite values only.
if ~all(isfinite([origin(:); ray(:); rlim(:); lambda(:); ...
        xgv(:); ygv(:); zgv(:)]))
    error('All input arguments must not be NaN or Inf.')
end

% Check whether RLIM has the correct number of elements.
if numel(rlim) ~= 2
    error('RLIM must have exactly 2 elements.')
end

% Check whether RLIM is ordered.
if diff(rlim) <= 0
    error('RLIM(2) must be greater than RLIM(1).');
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv)

% Check whether lambda has the correct size.
if any(size(lambda) ~= [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1])
    error('Size of LAMBDA does not match grid vectors.')
end

%% Compute probability of NaN measurements.
% Compute the ray from the origin to the maximum sensor range.
ray = ray / norm(ray) * rlim(2);

% Determine the number of rays.
nray = size(origin, 1);

% Preallocate the return matrix.
p = zeros(nray, 1);

% Compute the line parameter from origin to minimum sensor range.
tmin = rlim(1)/rlim(2);

% Loop over all rays and compute the respective probabilities for NaN
% measurements.
parfor i = 1 : nray
    % Compute the indices of the grid cells that the ray traverses from the
    % origin to the maximum sensor range.
    [vi, t] = trav(origin(i,:), ray(i,:), xgv, ygv, zgv);
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(size(lambda), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the length of the ray apportioned to each voxel when
    % traversing the grid from origin to maximum sensor range.
    d = diff(t) * norm(ray(i,:));
    
    % Compute the probability of obtaining an NaN measurement between 
    % maximum sensor range and infinity.
    psup = exp(-sum(lambda(vi) .* d));
    
    % Compute the length of the ray apportioned to each voxel when
    % traversing the grid from origin to minimum sensor range.
    t = [t(t < tmin); tmin];
    d = diff(t) * norm(ray(i,:));
    
    % Compute the probability of obtaining an NaN measurement between
    % origin and minimum sensor range.
    psub = 1 - exp(-sum(lambda(vi(1:max([length(d), 1]))) .* d));
    
    % Compute the overall probability of obtaining an NaN measurement.
    p(i) = psub + psup;
end

end
