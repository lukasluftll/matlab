function [p, L] = decayray(ls, lambda)
% DECAYRAY Compute probability of laser scan given ray decay map.
%   [P, L] = DECAYRAY(LS, LAMBDA) computes the probability of obtaining the 
%   laser scan LS conditioned on the decay map LAMBDA.
%
%   LS is a laserscan object containing N rays. The sensor pose of the scan
%   is assumed to be specified with respect to the decay rate map frame.
%
%   LAMBDA is a voxelmap object that contains the mean decay rate of each 
%   map voxel.
%
%   P and L are N-element row vectors. Together, they indicate the 
%   measurement probability for each ray.
%
%   If the m-th ray is a no-return, P(m) gives the corresponding 
%   measurement probability. For no-return rays, L(m) is zero.
%
%   If the m-th ray is reflected back to the sensor, the measurement 
%   probability is expressed by L(m). L(m) is the logarithm of the 
%   measurement probability density along the ray, evaluated at the ray 
%   endpoint. For returned rays, P(m) is unity.
%
%   Example:
%      ls = lsread('data/pcd/campus/pcd_sph/campus-00001.pcd', [2,120]);
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
% Compute the Cartesian ray direction with respect to the map frame.
ray = dir2cart(ls);

% Compute the logical indices of the returned rays.
iret = ls.ret;

% Set the length of no-return rays to maximum sensor range.
radius = ls.radius;
radius(~iret) = ls.rlim(2);
ray = ray .* repmat(radius, 1, 3);

%% Compute measurement probability for all rays.
% Preallocate the result matrices.
p = ones(ls.count, 1);
L = zeros(ls.count, 1);

% Compute the line parameter from origin to minimum sensor range.
tmin = ls.rlim(1) / ls.rlim(2);

% Loop over all rays.
parfor i = 1 : ls.count
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(tform2trvec(ls.sp(:,:,i)), ray(i,:), ...
        lambda.xgv, lambda.ygv, lambda.zgv); %#ok<PFBNS>

    % Convert the subscript voxel indices to linear indices.
    vi = sub2ind(size(lambda.data), vi(:,1), vi(:,2), vi(:,3));

    % Compute the length of the ray apportioned to each voxel.
    d = diff(t) * radius(i);
    
    % Compute log-likelihood or probability of the ray measurement.
    if iret(i) % Ray returned.
        % Compute the log-likelihood of the returned ray.
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
