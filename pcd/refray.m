function p = refray(origin, ray, rlim, ref)
% REFRAY Compute probability of Lidar measurement from reflectivity map.
%   P = REFRAY(ORIGIN, RAY, RLIM, REF) computes the probability of 
%   obtaining the Lidar measurement defined by ORIGIN and RAY conditioned 
%   on the reflectivity map REF.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the Cartesian 
%   origin coordinates and the ray vectors of the M measured rays. If all 
%   rays originate from the same point, ORIGIN may also be a 1x3 matrix.
%
%   RLIM is a 2-element vector that defines the minimum and the maximum
%   radius detected by the Lidar sensor. RAY rows whose lengths are not 
%   element of the interval defined by RLIM are assumed to be no-return 
%   measurements. For these measurements, RAY carries only information 
%   about the direction of the ray, not about its length.
%
%   REF is a voxelmap object that contains the reflectivity of each map 
%   voxel.
%
%   P is an M-element column vector. The value of the m-th element
%   corresponds to the probability of obtaining the m-th measurement.
%
%   Example:
%      origin = [0, 0, 0];
%      ray = [3, 4, 5];
%      rlim = [1, 10];
%      ref = voxelmap(repmat(magic(5)/100, [1, 1, 5]), 0:5, 0:5, 0:5);
%      p = refray(origin, ray, rlim, ref)
%
%   See also VOXELMAP, REFMAP, DECAYRAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(4, 4)

% Check if the arguments have the expected numbers of dimensions.
if ~ismatrix(origin) || ~ismatrix(ray)
    error('ORIGIN and RAY must be 2D matrices.')
end

% Check if the arguments have the expected sizes.
if size(origin, 2) ~= 3 || size(ray, 2) ~= 3
    error('ORIGIN and RAY must have 3 columns.')
end

% If ORIGIN has only one row, expand it to match the row size of RAY.
if size(origin, 1) == 1
    origin = repmat(origin, size(ray, 1), 1);
end

% Make sure all input arguments contain finite values only.
if ~all(isfinite([origin(:); ray(:); rlim(:); ref.data(:)]))
    error('Input arguments must not be NaN or Inf.')
end

% Check whether RLIM is ordered.
if diff(rlim) <= 0
    error('RLIM(2) must be greater than RLIM(1).');
end

%% Preprocess input arguments.
% Compute ray lengths.
l = sqrt(sum(ray.^2, 2));

% Determine which rays are no-returns.
inan = l < rlim(1) | l > rlim(2);

% Set the length of no-return rays to maximum sensor range.
if any(inan)
    ray(inan,:) = ray(inan,:) ./ repmat(l(inan), 1, 3) * rlim(2);
end

%% Compute probability of measurements.
% Loop over all rays.
nray = size(ray, 1);
p = zeros(nray, 1);
parfor i = 1 : nray
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(origin(i,:), ray(i,:), ...
        ref.xgv, ref.ygv, ref.zgv); %#ok<PFBNS>
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(size(ref.data), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the measurement probability depending on whether or not the 
    % ray returned.
    if inan(i) % Ray does not return.
        % Compute the indices of the voxels on the ray directly in front
        % and behind the measurement interval.
        ilim = knnsearch(t*rlim(2), rlim(:)) - 1; %#ok<PFBNS>
        
        % Calculate the probability that the ray is reflected before 
        % reaching the minimum sensor range.
        isub = vi(1 : ilim(1));
        psub = 1 - prod(1-ref.data(isub));
    
        % Calculate the probability that the ray surpasses the maximum 
        % sensor range.
        isup = vi(1 : ilim(2));
        psup = prod(1-ref.data(isup));
    
        % Sum up the probabilities to get the probability of the ray being
        % reflected before or after the measurement interval.
        p(i) = psub + psup;
    else % Ray returns.
        % Compute the probability of the ray being reflected in the last
        % voxel it traverses.
        p(i) = ref.data(vi(end)) * prod(1 - ref.data(vi(1:end-1)));
    end
end

end
