function [p, L] = decayray(ls, lambda)
% DECAYRAY Compute probability of laser scan given ray decay map.
%   [P, L] = DECAYRAY(LS, LAMBDA) computes the probability of obtaining the 
%   laser scan LS conditioned on the decay map LAMBDA.
%
%   LS is a laserscan object. The sensor pose of the scan is assumed to be 
%   specified with respect to the decay rate map frame.
%
%   LAMBDA is a voxelmap object that contains the mean decay rate of each 
%   map voxel.
%
%   The measurement probability is computed for each ray individually. 
%   If the m-th ray is reflected back to the sensor, the measurement 
%   probability is expressed by L(m). L(m) is the logarithm of the 
%   measurement probability density along the ray, evaluated at the ray 
%   endpoint. For returned rays, P(m) is unity.
%   If the m-th ray is a no-return, P(m) directly indicates the 
%   corresponding measurement probability. For no-return rays, L(m) is 
%   zero.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      ls = laserscan(pcd.azimuth, pcd.elevation, pcd.radius, [1, 100]);
%      lambda = decaymap(ls, -100:5:100, -100:5:100, -20:5:20);
%      [p, L] = decayray(ls, lambda)
%
%   See also LASERSCAN, VOXELMAP, DECAYMAP, REFRAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
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

% Compute the logical indices of the returned rays.
iret = ret(ls);

% Set the length of no-return rays to maximum sensor range.
ray(~iret,:) = ray(~iret,:) * ls.rlim(2);

%% Compute measurement probability for all rays.
% Preallocate the result matrices.
p = ones(ls.count, 1);
L = zeros(ls.count, 1);

% Compute the line parameter from origin to minimum sensor range.
tmin = ls.rlim(1) / ls.rlim(2);

% Loop over all rays.
parfor i = 1 : ls.count
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(ls.position, ray(i,:), ...
        lambda.xgv, lambda.ygv, lambda.zgv); %#ok<PFBNS>

    % Convert the subscript voxel indices to linear indices.
    vi = sub2ind(size(lambda.data), vi(:,1), vi(:,2), vi(:,3));

    % Compute the length of the ray apportioned to each voxel.
    d = diff(t) * norm(ray(i,:));
    
    % Compute log-likelihood or probability of ray measurement.
    if iret(i) % Ray returned.
        % Compute log-likelihood of returned rays.
        L(i) = log(lambda.data(vi(end))) - sum(lambda.data(vi) .* d);
    else % No-return ray.
        % Compute the probability of obtaining a no-return measurement
        % between maximum sensor range and infinity.
        psup = exp(-sum(lambda.data(vi) .* d));

        % Compute the length of the ray apportioned to each voxel when
        % traversing the grid from origin to minimum sensor range.
        t = [t(t < tmin); tmin];
        d = diff(t) * ls.rlim(2);

        % Compute the probability of obtaining an NaN measurement between
        % origin and minimum sensor range.
        psub = 1 - exp(-sum(lambda.data(vi(1:max([length(d), 1]))) .* d));

        % Compute the overall probability of obtaining an NaN measurement.
        p(i) = psub + psup;
    end
end

end
