function [lambda, r, l] = decaymap(ls, xgv, ygv, zgv)
% DECAYMAP Compute decay rate map from laser scans in grid volume.
%   LAMBDA = DECAYMAP(LS, XGV, YGV, ZGV) uses the laser scans LS to compute 
%   the mean ray decay rate LAMBDA for each voxel in the grid volume 
%   defined by the grid vectors XGV, YGV, ZGV.
%
%   LS is vector of laserscan objects. The sensor poses of the scans are
%   assumed to be specified with respect to the decay rate map frame.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   LAMBDA is a voxelmap object that contains the mean decay rate of each 
%   voxel. The lambda value of a voxel that has not been visited by any ray
%   is NaN.
%
%   [LAMBDA, R, L] = DECAYMAP(LS, XGV, YGV, ZGV) also returns the voxelmap 
%   objects R and L. 
%   R contains the number of ray remissions for each voxel. 
%   L contains the cumulated length of all rays that traversed the 
%   respective grid cell.
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
%      pcd = pcdread('castle.pcd');
%      ls = laserscan(pcd.azimuth, pcd.elevation, pcd.radius);
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
    warning(['LS.RLIM(1) > 0, but all no-return ray lengths are ', ...
      'assumed to surpass LS.RLIM(2), not to fall into [0; LS.RLIM(1)].'])
end

%% Compute ray lengths and returns per voxel.
% Determine the size of the voxel grid.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;

% Loop over all laser scans.
for s = 1 : numel(ls)
    % Compute the Cartesian ray direction vectors in the sensor frame.
    ray = dir2cart(ls(s));

    % Compute the indices of the returned rays.
    iret = ret(ls(s));

    % Set the length of no-return rays to maximum sensor range.
    ray(iret,:) = ray(iret,:) .* repmat(ls(s).radius(iret),1,size(ray,2));
    ray(~iret,:) = ray(~iret,:) * ls(s).rlim(2);

    % Use multiple workers to compute ray length and number of returns per 
    % voxel.
    spmd
        % Create return matrices for this worker.
        lw = zeros(gridsize);
        rw = zeros(gridsize);

        % For all rays of the worker's share compute the ray length per 
        % voxel.
        for i = labindex : numlabs : ls(s).count   
            % Compute the indices of the voxels through which the ray 
            % travels.
            pos = tform2trvec(ls(s).sp(:,:,i));
            [vi, t] = trav(pos, ray(i,:), xgv, ygv, zgv);

            % Convert the subscript indices to linear indices.
            vi = sub2ind(gridsize, vi(:,1), vi(:,2), vi(:,3));

            % Sum up the lengths the rays travel in each voxel.
            lw(vi) = lw(vi) + diff(t) * norm(ray(i,:));

            % In case of reflection, increment the number of returns.
            rw(vi(end)) = rw(vi(end)) + (iret(i) && t(end)==1);
        end
    end

    % Merge the results of all workers.
    r = voxelmap(rw{1}, xgv, ygv, zgv);
    l = voxelmap(lw{1}, xgv, ygv, zgv);
    for w = 2 : numel(lw)
        r = r + rw{w};
        l = l + lw{w};
    end
end

%% Compute the decay rate.
lambda = r ./ l;
lambda.data(l == 0) = NaN;

end
