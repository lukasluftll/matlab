function L = decayray(origin, ray, lambda)
% DECAYRAY Compute log-likelihood of Lidar measurement given ray decay map.
%   L = DECAYRAY(ORIGIN, RAY, LAMBDA, XGV, YGV, ZGV) computes the
%   log-likelihood of obtaining the Lidar ray measurement defined by ORIGIN 
%   and RAY conditioned on the decay map LAMBDA with grid vectors XGV, YGV,
%   ZGV.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the Cartesian 
%   origins and ray vectors of the M measured rays. If all rays originate
%   from the same point, ORIGIN may also be a 1x3 matrix.
%
%   LAMBDA is a voxelmap object that contains the mean decay rate of each 
%   map voxel.
%
%   L is an M-element column vector. The value of the m-th element
%   corresponds to the log-likelihood of obtaining the m-th measurement.
%   L is not equal to the log-likelihood of the observation, but it is
%   shifted by an unknown offset.
%
%   Example:
%      origin = [0, 0, 0];
%      ray = [3, 4, 5];
%      lambda = voxelmap(repmat(magic(5)/100, [1, 1, 5]), 0:5, 0:5, 0:5);
%      L = decayray(origin, ray, lambda)
%
%   See also VOXELMAP, DECAYNANRAY, DECAYMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(3, 3)

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
if ~all(isfinite([origin(:); ray(:); lambda.data(:)]))
    error('Input arguments must not be NaN or Inf.')
end

%% Compute log-likelihood of measurements.
% Determine the number of rays.
nray = size(origin, 1);

% Preallocate the return matrix.
L = zeros(nray, 1);

% Loop over all rays.
parfor i = 1 : nray
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(origin(i,:), ray(i,:), ...
        lambda.xgv, lambda.ygv, lambda.zgv); %#ok<PFBNS>
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(size(lambda.data), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the length of the ray apportioned to each voxel.
    d = diff(t) * norm(ray(i,:));

    % Compute the log-likelihood of the measurement.
    L(i) = log(lambda.data(vi(end))) - sum(lambda.data(vi) .* d);
end

end
