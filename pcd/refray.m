function p = refray(mts, ray, rlim, ref)
% REFRAY Compute probability of Lidar measurement from reflectivity map.
%   P = REFRAY(MTS, RAY, RLIM, REF) computes the probability of obtaining 
%   the Lidar measurement defined by ORIGIN and RAY conditioned on the 
%   reflectivity map REF.
%
%   MTS is an affine3d object that defines the pose of the sensor with
%   respect to the reflectivity map frame.
%
%   RAY is a Mx3 matrix whose rows contain the ray direction vectors of the 
%   M measured rays in the sensor frame.
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
%      ray = [3, 4, 5];
%      rlim = [1, 10];
%      ref = voxelmap(repmat(magic(5)/100, [1, 1, 5]), 0:5, 0:5, 0:5);
%      p = refray(affine3d(), ray, rlim, ref)
%
%   See also AFFINE3D, VOXELMAP, REFMAP, DECAYRAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(4, 4)

% Check the sensor pose.
if ~isa(mts, 'affine3d')
    error('MTS must be an affine3d object.')
end

% Check the ray matrix.
if ~ismatrix(ray)
    error('RAY must be a 2D matrix.')
end
if size(ray, 2) ~= 3
    error('RAY must have 3 columns.')
end

% Make sure all input arguments contain finite values only.
if ~all(isfinite([ray(:); rlim(:); ref.data(:)]))
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

% Set the length of no-return rays to maximum sensor range plus the
% diameter of the largest voxel.
rnan = rlim(2) + sqrt(3)*max([diff(ref.xgv),diff(ref.ygv),diff(ref.zgv)]);
ray(inan,:) = ray(inan,:) ./ [zeros(0, 3); repmat(l(inan), 1, 3)] * rnan;

% Transform rays from the sensor frame to the map frame.


%% Compute probability of measurements.
% Loop over all rays.
nray = size(ray, 1);
p = zeros(nray, 1);
parfor i = 1 : nray
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(mts(i,:), ray(i,:), ...
        ref.xgv, ref.ygv, ref.zgv); %#ok<PFBNS>
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(size(ref.data), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the measurement probability depending on whether or not the 
    % ray returned.
    if inan(i) % Ray does not return.
        % Compute the indices of the voxels on the ray directly in front
        % of the beginning and in front of the end of the measurement 
        % interval.
        ilim = knnsearch(t*rnan, rlim(:)) - 1; %#ok<PFBNS>
        
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
