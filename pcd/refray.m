function p = refray(ls, ref)
% REFRAY Compute probability of Lidar scan from reflectivity map.
%   P = REFRAY(LS, REF) computes the probability of obtaining the Lidar 
%   scan LS conditioned on the reflectivity map REF.
%
%   LS is a laserscan object. The sensor pose of the scan is assumed to be 
%   specified with respect to the reflectivity map frame.
%
%   REF is a voxelmap object that contains the reflectivity of each map 
%   voxel.
%
%   P is an M-element column vector. The value of the m-th element
%   corresponds to the probability of obtaining the m-th measurement.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      ls = laserscan(pcd.azimuth,pcd.elevation,pcd.radius,eye(4),[1,100]);
%      ref = refmap(ls, 0:5, 0:5, 0:5);
%      p = refray(ls, ref)
%
%   See also LASERSCAN, VOXELMAP, REFMAP, DECAYRAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(2, 2)

% Check the types of the input arguments.
if ~isa(ls, 'laserscan')
    error('LS must be a laserscan object.')
end
if ~isa(ref, 'voxelmap')
    error('REF must be a voxelmap object.')
end

%% Preprocess input arguments.
% Set the length of no-return rays to maximum sensor range plus the
% diameter of the largest voxel.
rnan = rlim(2) + sqrt(3)*max([diff(ref.xgv),diff(ref.ygv),diff(ref.zgv)]);
radius = ls.radius;
radius(noret(ls)) = rnan;

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
