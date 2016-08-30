function [lambda, r, l] = decaymap(ls, xgv, ygv, zgv)
% DECAYMAP Compute decay rate map from laser scans in grid volume.
%   LAMBDA = DECAYMAP(LS, XGV, YGV, ZGV) uses the laser scans LS to compute 
%   the mean ray decay rate LAMBDA for each voxel in the grid volume 
%   defined by the grid vectors XGV, YGV, ZGV.
%
%   LS is a vector of laserscan objects. The sensor poses of the scans are
%   assumed to be specified with respect to the decay rate map frame.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i,j,k] contains all points [x,y,z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   LAMBDA is a voxelmap object that contains the mean decay rate of each 
%   voxel. The lambda value of a voxel that has not been visited by any ray
%   is NaN. The prior of LAMBDA is the quotient of the number of returns 
%   by the total length of all rays.
%
%   [LAMBDA, R, L] = DECAYMAP(LS, XGV, YGV, ZGV) also returns the voxelmap 
%   objects R and L. 
%   R contains the number of ray remissions for each voxel. 
%   L contains the cumulated length of all rays that traversed the 
%   respective grid cell.
%   The priors of R and L are the total number of returns and the total 
%   length of all rays, respectively.
%
%   Concept of ray decay rate
%   -------------------------
%   The decay rate of a ray emitted by a Lidar sensor is a property of the
%   material through which the ray travels. The mean decay rate over a 
%   voxel is the number of ray returns from inside the voxel divided by the 
%   sum of the lengths of all rays travelling through the voxel:
%
%                  n_returns
%      lambda = -----------------
%                sum(ray_length)
%
%   The higher the sum of the ray lengths, the more accurate the
%   approximation of the decay rate.
%
%   Example:
%      ls = lsread('pcd/data/campus/pcd_sph/campus-00100.pcd');
%      lambda = decaymap(ls, -100:5:100, -100:5:100, -20:5:20)
%
%   See also LASERSCAN, VOXELMAP, DECAYRAY, REFMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(4, 4);

% Check the laser scan.
if ~isa(ls, 'laserscan')
    error('LS must be a laserscan object.')
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

% If the sensor measurement range starts at a positive value, issue a
% warning.
if ls.rlim(1) > 0
    warning('decaymap:rlim', ...
        ['LS.RLIM(1) > 0, but all no-return ray lengths are assumed ', ...
        'to surpass LS.RLIM(2), not to fall into [0; LS.RLIM(1)].'])
end

%% Compute ray lengths and returns per voxel.
% Determine the size of the voxel grid.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;

% Loop over all laser scans.
ltot = 0;
rtot = 0;
for s = 1 : numel(ls)
    % Compute the Cartesian ray direction vectors in the map frame.
    ray = dir2cart(ls(s));

    % Set the length of no-return rays to maximum sensor range.
    radius = ls(s).radius;
    radius(~ls(s).ret) = ls(s).rlim(2);
    ray = ray .* repmat(radius, 1, 3);
    
    % Sum up the lengths of all rays and the number of returns.
    ltot = ltot + sum(radius);
    rtot = rtot + sum(ls(s).ret);

    % Compute ray length and number of returns per voxel.
    l = zeros(gridsize);
    r = zeros(gridsize);
    for i = 1 : ls(s).count   
        % Compute the indices of the voxels through which the ray travels.
        [vi,t] = trav(tform2trvec(ls(s).sp(:,:,i)), ray(i,:), xgv,ygv,zgv);

        % Convert the subscript indices to linear indices.
        vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));

        % Sum up the lengths the rays travel in each voxel.
        l(vi) = l(vi) + diff(t) * radius(i);

        % In case of reflection, increment the number of returns.
        r(vi(end)) = r(vi(end)) + (iret(i) && t(end)==1);
    end
end

%% Create voxelmaps.
r = voxelmap(r, xgv, ygv, zgv, rtot);
l = voxelmap(l, xgv, ygv, zgv, ltot);
lambda = r ./ l;
lambda.data(l==0) = NaN;

end
