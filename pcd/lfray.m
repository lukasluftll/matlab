function [p, L] = lfray(ls, lf, pnr)
% LFRAY Compute probability of laser scan given likelihood field.
%   [P, L] = LFRAY(LS, LF, PNR) computes the probability of obtaining the 
%   laser scan LS conditioned on the likelihood field map LF.
%
%   LS is a laserscan object that contains N rays. The sensor pose of the
%   scan is assumed to be specified with respect to the likelihood field
%   map frame.
%
%   LF is a voxelmap object that contains the likelihood field. The
%   likelihood field is a measure for the probability that a map voxel
%   belongs to an object. It may or may not be normalized; LFRAY normalizes 
%   the return probability over each ray.
%
%   PNR defines the unconditioned probability of obtaining a no-return ray.
%   Usually, it is defined as the ratio of the number of no-return rays to 
%   the number of all rays.
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
%   If the points on a returned ray corresponding to minimum and maximum 
%   sensor range lie outside the map, the probabilities cannot be 
%   normalized. The respective values of p and L are set to NaN.
%
%   Example:
%      ls = lsread('pcd/data/sph.pcd', [2,120]);
%      lf = lfmap(ls2pc(ls), 1, -100:5:100, -100:5:100, -20:5:20);
%      [p, L] = lfray(ls, constrain(lf, [0.001,Inf]), 0.2)
%
% See also LASERSCAN, VOXELMAP, LFMAP, LFPC, DECAYRAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(3, 3)

% Check input arguments.
validateattributes(ls, {'laserscan'}, {}, '', 'LS')
validateattributes(lf, {'voxelmap'}, {}, '', 'LF')
validateattributes(pnr, {'numeric'}, {'>=',0, '<=',1}, '', 'PNR')

%% Compute scan probability.
% Compute the logical indices of the returned rays.
iret = ls.ret;

% Loop over all rays and compute their probabilities.
p = ones(ls.count, 1);
L = zeros(ls.count, 1);
parfor i = 1 : ls.count
    if iret(i) % Returned ray.
        % Compute the indices of the grid cells that the ray traverses
        % inside the valid sensor range.
        ray = dir2cart(ls, i); %#ok<*PFBNS>
        start = tform2trvec(ls.sp(:,:,i)) + ray*ls.rlim(1);
        [iv,t] = trav(start, ray*diff(ls.rlim), lf.xgv, lf.ygv, lf.zgv); 
        
        % If the map does not cover the segment of the ray from minimum to
        % maximum sensor range, set the probability of the ray to NaN.
        if isempty(t) || t(1) > 0 || t(end) < 1
            p(i) = NaN;
            L(i) = NaN;
            continue
        end
        
        % Convert the subscript voxel indices to linear indices.
        iv = sub2ind(size(lf.data), iv(:,1), iv(:,2), iv(:,3));
        
        % Compute the non-normalized probability of a return.
        pr = sum(lf.data(iv) .* diff(t) * diff(ls.rlim));

        % Compute the normalized log-likelihood of obtaining the returned 
        % ray.
        L(i) = log(lf.data(iv(end))) + log(1-pnr) - log(pr);
    else % No-return ray.
        p(i) = pnr;
    end
end

end
