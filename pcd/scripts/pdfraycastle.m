% Visualizes the probability of obtaining a Lidar scan that is shifted with
% respect to itself.
%
% Steps:
% # Loads the castle scan.
% # Computes a Lidar decay map from it.
% # Shifts the scan horizontally in x and y direction.
% # Computes the probability of obtaining the shifted scan with respect to 
%   the decay map.

% Read the point cloud.
pcd = pcdread('data/castle.pcd');

% Determine the number of rays (both returns and no-returns).
nray = numel(pcd.azimuth);

% Define the grid vectors for voxelization.
res = 1;
xgv = min(pcd.x(:)) : res : max(pcd.x(:));
ygv = min(pcd.y(:)) : res : max(pcd.y(:));
zgv = min(pcd.z(:)) : res : max(pcd.z(:));

% Determine the number of returns.
npoints = sum(isfinite(pcd.radius(:)));

% Compute the decay rate prior.
maxrange = 100;   % VLP-16 can see 100 m far.
raylength = pcd.radius;
raylength(~isfinite(raylength)) = maxrange;
lambdaPrior = npoints / sum(raylength(:));

% Compute the ray direction vectors.
[dirx, diry, dirz] = sph2cart(...
    pcd.azimuth(:), pcd.elevation(:), pcd.radius(:));

% Compute the decay rate map.
lambdaLim = [5e-4, 1e1];
lambda = raydecay(pcd.azimuth, pcd.elevation, pcd.radius, ...
    xgv, ygv, zgv);

% Limit the decay rates.
lambda = max(lambdaLim(1), lambda);
lambda = min(lambdaLim(2), lambda);

% Add assumption about decay rate prior.
lambda(~isfinite(lambda)) = lambdaPrior;

% Shift the scan and compute the probability of obtaining it.
xs = -3 : 1 : 3;
ys = -3 : 1 : 3;
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

% Display the overall probabilities of the shifted scans.
surf(xs, ys, prob);
