% Loads the castle scan, computes a Lidar decay map from it, and then 
% shifts the scan horizontally and computes the probability of obtaining 
% the shifted scan with respect to the decay map.

% Read the point cloud.
pcd = pcdread('data/castle.pcd');

% Determine the number of rays.
nray = numel(pcd.azimuth);

% Define the grid vectors for voxelization.
res = 5;
xgv = min(pcd.x(:)) : res : max(pcd.x(:));
ygv = min(pcd.y(:)) : res : max(pcd.y(:));
zgv = min(pcd.z(:)) : res : max(pcd.z(:));

% Replace no-return ray lengths with a value larger than the grid diameter.
pcd.radius(~isfinite(pcd.radius)) = 1e4;

% Compute the ray direction vectors.
[dirx, diry, dirz] = sph2cart(...
    pcd.azimuth(:), pcd.elevation(:), pcd.radius(:));

% Compute the decay rate map.
lambda = raydecay(pcd.azimuth, pcd.elevation, pcd.radius, ...
    xgv, ygv, zgv);
lambda(isnan(lambda)) = 0;

% Shift the scan and compute the probability of obtaining it.
xs = -0.1 : 0.1 : 0.1;
ys = -0.1 : 0.1 : 0.1;
prob = zeros(length(xs), length(ys));
waitbarHandle = waitbar(0, 'Computing scan probabilities ...');
for i = 1 : length(xs)
    for j = 1 : length(ys)
        origin = [xs(i), ys(j), 0];
        p = pdfray(origin, [dirx, diry, dirz], lambda, xgv, ygv, zgv);
        prob(i,j) = sum(log(p));
        
        % Advance the progress bar.
        waitbar(((i-1)*length(xs) + j) / (length(xs)*length(ys)), ...
            waitbarHandle);
    end
end
close(waitbarHandle);
