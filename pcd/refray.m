function [p, L] = refray(ls, ref)
% REFRAY Compute probability of Lidar scan from reflectivity map.
%   P = REFRAY(LS, REF) computes the probability of obtaining the Lidar 
%   scan LS conditioned on the reflectivity map REF.
%
%   LS is a laserscan object. The sensor pose of the scan is assumed to be 
%   specified with respect to the reflectivity map coordinate frame.
%
%   REF is a voxelmap object that contains the reflectivity of each map 
%   voxel.
%
%   P is an M-element column vector. The value of the m-th element
%   corresponds to the probability of obtaining the m-th measurement.
%
%   [P, L] = REFRAY(LS, REF) also returns the M-element column vector L.
%   L(m) is the logarithm of the probability density of the m-th ray, 
%   evaluated at the ray endpoint. For no-return rays, L(m) is zero.
%
%   Example:
%      pcd = pcdread('castle.pcd');
%      ls = laserscan(pcd.azimuth, pcd.elevation, pcd.radius, [1,100]);
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
% Compute the Cartesian ray direction vectors.
ray = cart(ls);

% Compute the indices of the returned rays.
iret = ret(ls);

% Set the length of no-return rays to maximum sensor range plus the
% diameter of the largest voxel.
radiusnr = ls.rlim(2) + ...
    sqrt(3) * max([diff(ref.xgv), diff(ref.ygv), diff(ref.zgv)]);
ray(~iret,:) = ray(~iret,:) * radiusnr;

%% Compute probability of measurements.
% Loop over all rays.
p = zeros(ls.count, 1);
L = zeros(ls.count, 1);
parfor i = 1 : ls.count
    % Compute the indices of the grid cells that the ray traverses.
    [vi,t] = trav(ls.position,ray(i,:),ref.xgv,ref.ygv,ref.zgv);%#ok<PFBNS>
    
    % Convert the subscript voxel indices to linear indices.
    vi = sub2ind(size(ref.data), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the measurement probability depending on whether or not the 
    % ray returned.
    if iret(i) % Ray returns.
        % Compute the probability of the ray being reflected in the last
        % voxel it traverses.
        p(i) = ref.data(vi(end)) * prod(1 - ref.data(vi(1:end-1)));
        
        % Compute the log-likelihood normalized over the ray length.
        % TODO: Make sure the ray traverses the voxel.
        L(i) = log(p(i) / (t(end) * norm(ray(i,:))));
    else % Ray does not return.
        % Compute the indices of the voxels on the ray directly in front
        % of the beginning and in front of the end of the measurement 
        % interval.
        ilim = knnsearch(t*radiusnr, ls.rlim(:)) - 1;
        
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
    end
end

end
