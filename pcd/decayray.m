function L = decayray(ls, lambda)
% DECAYRAY Compute log-likelihood of laser scan given ray decay map.
%   L = DECAYRAY(LS, LAMBDA) computes the log-likelihood of obtaining the 
%   laser scan LS conditioned on the decay map LAMBDA.
%
%   LS is a laserscan object. The sensor pose of the scan is assumed to be 
%   specified with respect to the decay rate map frame.
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
%      pcd = pcdread('castle.pcd');
%      ls = laserscan(pcd.azimuth, pcd.elevation, pcd.radius, [1, 100]);
%      lambda = decaymap(ls, -100:5:100, -100:5:100, -20:5:20);
%      L = decayray(ls, lambda)
%
%   See also LASERSCAN, VOXELMAP, DECAYNANRAY, DECAYMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(2, 2)

% Check the input argument types.
if ~isa(ls, 'laserscan')
    error('LS must be a laserscan object.')
end
if ~isa(lambda, 'voxelmap')
    error('LAMBDA must be a voxelmap object.')
end

%% Preprocess input arguments.
% Compute the Cartesian ray direction vectors.
ray = cart(ls);

% Compute the indices of the returned rays.
iret = ret(ls);

% Set the length of no-return rays to maximum sensor range plus the
% diameter of the largest voxel.
radiusnr = ls.rlim(2) + ...
    sqrt(3) * max([diff(ref.xgv), diff(ref.ygv), diff(ref.zgv)]);
ray(~iret,:) = ray(~iret,:) * radiusnr;

%% Compute log-likelihood of measurements.
% Preallocate the return matrix.
L = zeros(ls.count, 1);

% Loop over all rays.
parfor i = 1 : ls.count
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(ls.position, ray(i,:), ...
        lambda.xgv, lambda.ygv, lambda.zgv); %#ok<PFBNS>
    
    % Convert the subscript voxel indices to linear indices.
    vi = sub2ind(size(lambda.data), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the length of the ray apportioned to each voxel.
    d = diff(t) * norm(ray(i,:));

    % Compute the log-likelihood of the measurement.
    L(i) = log(lambda.data(vi(end))) - sum(lambda.data(vi) .* d);
end

end
