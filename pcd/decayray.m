function [L, p] = decayray(ls, lambda)
% DECAYRAY Compute probability of laser scan given ray decay map.
%   [L, p] = DECAYRAY(LS, LAMBDA) computes the probability of obtaining the 
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
%   probability is expressed as log-likelihood L(m). L(m) is not 
%   equal to the log-likelihood of the ray, but it is shifted by an unknown 
%   offset. p(m) is unity.
%   If the m-th ray is a no-return, p(m) contains the corresponding 
%   measurement probability. L(m) is zero.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      ls = laserscan(pcd.azimuth, pcd.elevation, pcd.radius, [1, 100]);
%      lambda = decaymap(ls, -100:5:100, -100:5:100, -20:5:20);
%      [L, p] = decayray(ls, lambda)
%
%   See also LASERSCAN, VOXELMAP, DECAYMAP.

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

% Compute the indices of the returned rays.
iret = ret(ls);

% Set the length of no-return rays to maximum sensor range.
ray(~iret,:) = ray(~iret,:) * ls.rlim(2);

%% Compute log-likelihood of returned rays.
% Loop over all returned rays.
retray = ray(:,iret);
Lr = zeros(sum(ret), 1);
parfor i = 1 : size(retray, 1)
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(ls.position, retray(i,:), ...
        lambda.xgv, lambda.ygv, lambda.zgv); %#ok<PFBNS>
    
    % Convert the subscript voxel indices to linear indices.
    vi = sub2ind(size(lambda.data), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the length of the ray apportioned to each voxel.
    d = diff(t) * norm(ray(i,:));

    % Compute the log-likelihood of the measurement.
    Lr(i) = log(lambda.data(vi(end))) - sum(lambda.data(vi) .* d);
end

% Construct the log-likelihood vector for all rays.
L = zeros(ls.count, 1);
L(iret) = Lr;

%% Compute probability of no-return rays.
% Preallocate the return matrix.
pnr = zeros(ls.count-sum(ret), 1);

% Compute the line parameter from origin to minimum sensor range.
tmin = rlim(1) / rlim(2);

% Loop over all rays and compute the respective probabilities for NaN
% measurements.
parfor i = 1 : size(pnr, 1)
    % Compute the indices of the grid cells that the ray traverses from the
    % origin to the maximum sensor range.
    [vi, t] = trav(ls.position, ray(i,:), ...
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
