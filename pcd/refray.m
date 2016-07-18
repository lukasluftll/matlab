function p = refray(origin, ray, rlim, ref, xgv, ygv, zgv)
% REFRAY Compute probability of Lidar measurement from reflectivity map.
%   L = REFRAY(ORIGIN, RAY, REF, XGV, YGV, ZGV) computes the probability
%   of obtaining the Lidar ray measurement defined by ORIGIN and RAY 
%   conditioned on the reflectivity map REF with grid vectors XGV, YGV,
%   ZGV.
%
%   ORIGIN and RAY are Mx3 matrices whose rows contain the Cartesian 
%   origins and ray vectors of the M measured rays. If all rays originate
%   from the same point, ORIGIN may also be a 1x3 matrix.
%
%   RLIM is a 2-element vector that defines the minimum and the maximum
%   radius detected by the Lidar sensor. RAY values that are not element of
%   the interval defined by RLIM are assumed to be no-return measurements.
%   For these measurements, RAY carries only information about the
%   direction of the ray, not about its length.
%
%   XGV, YGV, ZGV are vectors that define the rasterization of the grid.
%   A voxel with index [i, j, k] contains all points [x, y, z] that satisfy
%   the inequality:
%
%      (XGV(i) <= x < XGV(i+1))
%      && (YGV(j) <= y < YGV(j+1)) 
%      && (ZGV(k) <= z < ZGV(k+1))
%
%   REF is a IxJxK matrix that contains the reflectivity of each map voxel,
%   where I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
%
%   L is an M-element column vector. The value of the m-th element
%   corresponds to the probability of obtaining the m-th measurement.
%
%   Example:
%      origin = [0, 0, 0];
%      ray = [3, 4, 5];
%      ref = repmat(magic(5)/100, [1, 1, 5]);
%      gv = 0 : 5; 
%      p = refray(origin, ray, lim, ref, gv, gv, gv)
%
%   See also REFMAP, DECAYRAY.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check whether the user provided the correct number of input arguments.
narginchk(7, 7)

% Check if the arguments have the expected numbers of dimensions.
if ~ismatrix(origin) || ~ismatrix(ray) || ndims(ref) ~= 3
    error('ORIGIN and RAY must be 2D matrices, REF must be 3D.')
end

% Check if the arguments have the expected sizes.
if size(origin, 2) ~= 3 || size(ray, 2) ~= 3
    error('ORIGIN and RAY must have 3 columns.')
end

% If ORIGIN has only one row, expand it to match the row size of RAY.
if size(origin, 1) == 1
    origin = repmat(origin, size(ray, 1), 1);
end

% Make sure all input arguments contain finite values only.
if ~all(isfinite([origin(:); ray(:); rlim(:); ref(:)]))
    error('Input arguments must not be NaN or Inf.')
end

% Check whether RLIM is ordered.
if diff(rlim) <= 0
    error('RLIM(2) must be greater than RLIM(1).');
end

% Check the grid vectors.
gvchk(xgv, ygv, zgv)

% Check whether REF has the correct size.
if any(size(ref) ~= [numel(xgv)-1, numel(ygv)-1, numel(zgv)-1])
    error('Size of REF does not match grid vectors.')
end

%% Compute probability of measurements.
% Determine the number of rays.
nray = size(origin, 1);

% Preallocate the return matrix.
p = zeros(nray, 1);

l = norm(ray);
inan = find(l < rlim(1) | l > rlim(2));

% Loop over all rays.
parfor i = 1 : nray
    % Compute the indices of the grid cells that the ray traverses.
    [vi, t] = trav(origin(i,:), ray(i,:), xgv, ygv, zgv);
    
    % Convert the subscript indices to linear indices.
    vi = sub2ind(size(ref), vi(:,1), vi(:,2), vi(:,3));
    
    % Compute the indices of the voxels on the ray that are next to the 
    % minimum and maximum measurement range barrier.
    [~, ilim] = min(abs(...
        repmat(t*norm(ray), 1, numel(rlim)) - repmat(rlim, numel(t), 1)));
    ilim = ilim - 1;
    
    % Calculate the probability that the ray is reflected before reaching
    % the minimum sensor range.
    isub = vi(1 : ilim(1));
    msub = ref(isub);
    psub = 1 - prod(1-msub);
    
    % Calculate the probability that the ray surpasses the maximum sensor 
    % range.
    isup = vi(ilim(1)+1 : ilim(2));
    msup = ref(isup);
    psup = prod(1-msup);
    
    % Sum up the probabilities.
    p(i) = psub + psup;
end

end
