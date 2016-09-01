function [p, L] = refray(ls, ref)
% REFRAY Compute probability of laser scan from reflectivity map.
%   [P, L] = REFRAY(LS, REF) computes the probability of obtaining the 
%   laser scan LS conditioned on the reflectivity map REF.
%
%   LS is a laserscan object containing N rays. The sensor pose of the scan
%   is assumed to be specified with respect to the reflectivity map 
%   coordinate frame.
%
%   REF is a voxelmap object that contains the reflectivity of each map 
%   voxel.
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
%   p(i) and L(i) are set to NaN if any of the following applies to ray i:
%   - The starting point lies outside the map.
%   - Ray i is a returned ray and its endpoint ray lies outside the map.
%   - Ray i is a no-return ray and the point on the ray corresponding to 
%     maximum sensor range lies outside the map.
%
%   Example:
%      ls = lsread('pcd/data/sph.pcd', [2,120]);
%      ref = refmap(ls, -100:5:100, -100:5:100, -20:5:20);
%      [p, L] = refray(ls, ref)
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

% Check the validity of the reflectivity map.
if any(ref.data(:) < 0) || any(ref.data(:) > 1)
    error('REF must contain values in [0;1].')
end

%% Preprocess input arguments.
% Compute the normalized Cartesian ray direction vectors with respect to
% the map frame.
ray = dir2cart(ls);

% Compute the logical indices of the returned rays.
iret = ret(ls);

% Set the length of no-return rays to maximum sensor range.
l = ls.radius;
l(~iret) = ls.rlim(2);

% Increase the length of each returned ray by the diameter of the largest 
% voxel to make sure it completely traverses the voxel that contains the
% ray endpoint.
l = l + sqrt(3)*max([diff(ref.xgv(:));diff(ref.ygv(:));diff(ref.zgv(:))]);
ray = ray .* repmat(l, 1, 3);

%% Compute probability of measurements.
% Loop over all rays.
p = ones(ls.count, 1);
L = zeros(ls.count, 1);
parfor i = 1 : ls.count
    % Compute the indices of the grid cells that the ray traverses.
    [iv,t] = trav(tform2trvec(ls.sp(:,:,i)), ray(i,:), ...
        ref.xgv, ref.ygv, ref.zgv); %#ok<PFBNS>
    
    % Convert the subscript voxel indices to linear indices.
    iv = sub2ind(size(ref.data), iv(:,1), iv(:,2), iv(:,3));
    
    % Compute the index of the voxel where the ray is reflected and the
    % index of the voxel where the sensor measurement range ends.
    ivret = find(t*l(i) < ls.radius(i), 1, 'last');
    ivend = find(t*l(i) < ls.rlim(2), 1, 'last');
    
    % If the starting point or the endpoint of the ray lie outside the map,
    % return NaN probability.
    if t(1) > 0 || (iret(i) && ivret >= numel(t)) ...
            || (~iret(i) && ivend >= numel(t))
        p(i) = NaN;
        L(i) = NaN;
        continue
    end
    
    % Compute the measurement probability depending on whether or not the 
    % ray returned.
    if iret(i) % Ray returns.      
        % Compute the length of the ray apportioned to this voxel.
        lr = diff(t(ivret + [0,1])) * l(i);
        
        % Compute the density of the probability that the ray is reflected.
        L(i) = log(ref.data(iv(ivret))) ...
            + sum(log(1-ref.data(iv(1:ivret-1)))) - log(lr);
    else % Ray does not return.
        % Compute the index of the voxel where the sensor measurement range
        % begins.
        ivstart = find(t*l(i) < ls.rlim(1), 1, 'last');
        
        % Compute the lengths of the ray apportioned to the voxels where 
        % the sensor measurement range starts and ends.
        dtstart = diff(t([ivstart; ivstart+1]));
        dtend = diff(t([ivend; ivend+1]));
        
        % Compute the weights for the reflectivity cells before the ray 
        % reaches minimum sensor range. The weights are chosen according to 
        % the ray lengths apportioned to the cells.
        w = [ones(1,ivstart-1), (ls.rlim(1)/l(i)-t(ivstart)) / dtstart];
        
        % Calculate the probability that the ray is reflected before 
        % reaching the minimum sensor range.
        psub = 1 - prod(1 - ref.data(1:ivstart) .* w);
        
        % Compute the weights for the reflectivity cells before the ray
        % reaches maximum sensor range.
        w = [ones(1, ivend-1), ...
            (min([ls.rlim(2)/l(i), t(end)]) - t(ivend)) / dtend];

        % Calculate the probability that the ray surpasses the maximum 
        % sensor range.
        psup = prod(1 - ref.data(1:ivend).*w);
    
        % Sum up the probabilities to get the probability of the ray being
        % reflected before or after the measurement interval.
        p(i) = psub + psup;
    end
end

end
