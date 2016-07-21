function [ix, iy, iz] = grdidx(p, xgv, ygv, zgv)
% GRDIDX Compute voxel index of 3D point.
%   [IX, IY, IZ] = GRDIDX(P, XGV, YGV, ZGV) computes the voxel indices IX, 
%   IY, IZ of the 3D points P in a voxel grid whose rasterization is 
%   defined by the grid vectors XGV, YGV, ZGV.
%
%   P is a Nx3 matrix whose rows contain the coordinates of the points.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the 
%   voxel grid. A voxel with index [i,j,k] contains all points [x,y,z]
%   that satisfy the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   IX, IY, IZ are N-element column vectors that contain the indices of the
%   voxels for every point. If a point does not reside inside the voxel
%   grid, it is assigned an index of 0 in IX, IY, and IZ.
%
%   Example:
%      [ix, iy, iz] = grididx([4, 5, 6], 1:10, 1:10, 1:10)

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(4, 4);

% Check the size of P.
if size(p, 2) ~= 3
    error('P must be a Nx3 matrix.')
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv);

%% Compute voxel indices.
% Use multiple workers to compute the voxel indices of all points.
spmd
    % Initialize the index matrix for this worker.
    iw = zeros(size(p));
    
    % Loop over this worker's share of all points.
    for i = labindex : numlabs : size(p,1)
        % Compute the index of the point.
        iwtmp = [find(xgv <= p(i,1), 1, 'last') * (p(i,1) < xgv(end)), ...
            find(ygv <= p(i,2), 1, 'last') * (p(i,2) < ygv(end)), ...
            find(zgv <= p(i,3), 1, 'last') * (p(i,3) < zgv(end))];

        % If the index could not be found, set it to 0.
        iw(i, [true,true,true] & numel(iwtmp)==3) = iwtmp;
    end
end

% Merge the results of all the workers.
idx = iw{1};
for i = 2 : numel(iw)
    idx = idx + iw{i};
end

ix = idx(:,1);
iy = idx(:,2);
iz = idx(:,3);

end
