% Read the point cloud.
pcd = pcdread('data/castle.pcd');

% Determine the number of rays.
nrays = numel(pcd.azimuth);

% Define the grid vectors for voxelization.
res = 5;
xgv = min(pcd.x(:)) : res : max(pcd.x(:));
ygv = min(pcd.y(:)) : res : max(pcd.y(:));
zgv = min(pcd.z(:)) : res : max(pcd.z(:));

% Replace no-return ray lengths with a value larger than the grid diameter.
pcd.radius(~isfinite(pcd.radius)) ...
    = sum(diff(xgv)) + sum(diff(ygv)) + sum(diff(zgv));

% Compute the ray direction vectors.
[dirx, diry, dirz] = sph2cart(pcd.azimuth(:), pcd.elevation(:), pcd.radius(:));

% Compute the decay rate map.
lambda = raydecay(pcd.azimuth, pcd.elevation, pcd.radius, ...
    xgv, ygv, zgv);

% Compute the probability of the scan.
origin = zeros(nrays, 3);
p = pdfray(origin, [dirx, diry, dirz], lambda, xgv, ygv, zgv);
