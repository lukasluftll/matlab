function p = decaynanray(origin, ray, rlim, lambda)
% DECAYNANRAY Compute probability of NaN Lidar measurement from decay map.
%   P = DECAYNANRAY(ORIGIN, RAY, RLIM, LAMBDA) computes the probability of 
%   obtaining the measurement NaN from a Lidar sensor that sends a ray from 
%   ORIGIN in direction RAY through a ray decay voxel grid defined by 
%   LAMBDA with grid vectors XGV, YGV, ZGV.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the Cartesian 
%   origins and ray vectors of the M measured rays. If all rays originate
%   from the same point, ORIGIN may also be a 1x3 matrix. The lengths of 
%   the rays are not considered.
%
%   RLIM is a 2-element vector that defines the minimum and the maximum
%   radius detected by the Lidar sensor. All other radii are assumed to
%   result in an NaN measurement.
%
%   LAMBDA is a voxelmap object that contains the mean decay rate of each 
%   map voxel.
%
%   P is an M-element column vector. The value of the m-th element
%   corresponds to the probability of obtaining NaN for the m-th 
%   measurement.
%
%   Example:
%      origin = [0, 0, 0];
%      ray = [3, 4, 5];
%      rlim = [1; 100];
%      lambda = voxelmap(repmat(magic(5)/100, [1, 1, 5]), 0:5, 0:5, 0:5);
%      p = decaynanray(origin, ray, rlim, lambda)
%
%   See also VOXELMAP, DECAYRAY, DECAYMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(4, 4);

% Check if the arguments have the expected sizes.
if ~(ismatrix(origin) && ismatrix(ray))
    error('ORIGIN and RAY must be 2D matrices.')
end
if size(origin, 2) ~= 3 || size(ray, 2) ~= 3
    error('ORIGIN and RAY must have 3 columns.')
end

% If ORIGIN has only one row, expand it to match the row size of RAY.
if size(origin, 1) == 1
    origin = repmat(origin, size(ray, 1), 1);
end

% Make sure all input arguments contain finite values only.
if ~all(isfinite([origin(:); ray(:); rlim(:); lambda.data(:)]))
    error('All input arguments must not be NaN or Inf.')
end

% Check whether RLIM has the correct number of elements.
if numel(rlim) ~= 2
    error('RLIM must have exactly 2 elements.')
end

% Check whether RLIM is sorted.
if diff(rlim) <= 0
    error('RLIM(2) must be greater than RLIM(1).');
end

%% Compute probability of NaN measurements.
% Compute the ray vector pointing from the origin to the maximum sensor 
% range.
ray = ray / norm(ray) * rlim(2);

% Determine the number of rays.
nray = size(origin, 1);

% Preallocate the return matrix.
p = zeros(nray, 1);

% Compute the line parameter from origin to minimum sensor range.
tmin = rlim(1) / rlim(2);

% Loop over all rays and compute the respective probabilities for NaN
% measurements.
parfor i = 1 : nray
    % Compute the indices of the grid cells that the ray traverses from the
    % origin to the maximum sensor range.
    [vi, t] = trav(origin(i,:), ray(i,:), ...
        lambda.xgv, lambda.ygv, lambda.zgv); %#ok<PFBNS>
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(size(lambda.data), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the length of the ray apportioned to each voxel when
    % traversing the grid from origin to maximum sensor range.
    d = diff(t) * rlim(2); %#ok<PFBNS>
    
    % Compute the probability of obtaining an NaN measurement between 
    % maximum sensor range and infinity.
    psup = exp(-sum(lambda.data(vi) .* d));
    
    % Compute the length of the ray apportioned to each voxel when
    % traversing the grid from origin to minimum sensor range.
    t = [t(t < tmin); tmin];
    d = diff(t) * rlim(2);
    
    % Compute the probability of obtaining an NaN measurement between
    % origin and minimum sensor range.
    psub = 1 - exp(-sum(lambda.data(vi(1:max([length(d), 1]))) .* d));
    
    % Compute the overall probability of obtaining an NaN measurement.
    p(i) = psub + psup;
end

end
