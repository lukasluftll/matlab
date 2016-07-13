% Visualizes the log-likelihood of obtaining a shifted Lidar scan given a
% map created from the same Lidar scan at the true pose.
%
% Steps:
% # Loads the scan.
% # Computes a Lidar decay map from it.
% # Shifts the scan horizontally in x and y direction.
% # Computes the log-likelihood of obtaining the shifted scan with respect
%   to the decay map.

%% Set parameters.
% Shifting offset in x and y direction.
shift = 3;

% Resolution of the decay rate map.
res = 1;

% Resolution of the log-likelihood graph.
shiftres = 1;

% Minimum and maximum range of the Lidar sensor.
rlim = [0, 130];

% Minimum and maximum admissible decay rate.
lambdaLim = [2e-3, 1e+1];

%% Compute decay rate map.
% Read the point cloud.
pcd = pcdread('data/castle.pcd');

% Compute the decay rate map.
radiusFinite = pcd.radius;
radiusFinite(~isfinite(radiusFinite)) = rlim(2);
xgv = min(pcd.x(:))-shift : res : max(pcd.x(:))+res+shift;
ygv = min(pcd.y(:))-shift : res : max(pcd.y(:))+res+shift;
zgv = min(pcd.z(:))-shift : res : max(pcd.z(:))+res+shift;
lambda = raydecay(pcd.azimuth, pcd.elevation, radiusFinite, xgv, ygv, zgv);

% Set all voxels without data to the decay rate prior.
lambda(~isfinite(lambda)) = ...
    sum(isfinite(pcd.radius(:))) / sum(radiusFinite(:));

% Limit the decay rates to a reasonable interval and add prior.
lambda = max(lambdaLim(1), lambda);
lambda = min(lambdaLim(2), lambda);

%% Compute log-likelihood of shifted scans.
% Compute the direction vectors of the returned rays.
[dirxr, diryr, dirzr] = sph2cart(pcd.azimuth(isfinite(pcd.radius)), ...
    pcd.elevation(isfinite(pcd.radius)), pcd.radius(isfinite(pcd.radius)));

% Compute the direction vectors of the no-return rays.
[dirxnr, dirynr, dirznr] = sph2cart(pcd.azimuth(~isfinite(pcd.radius)), ...
    pcd.elevation(~isfinite(pcd.radius)), 1);

% Shift the scan and compute the probability of obtaining it.
gvs = -shift : shiftres : shift;
prob = zeros(numel(gvs));
waitbarHandle = waitbar(0, 'Computing scan probabilities ...');
L = zeros(numel(gvs));
for i = 1 : numel(gvs)
    for j = 1 : numel(gvs)
        origin = [gvs(i), gvs(j), 0];
        
        % Compute the log-likelihood of the measurements, depending on
        % whether or not the individual ray returned.
        Lr = sum(pdfray(origin, [dirxr, diryr, dirzr], lambda, ...
            xgv, ygv, zgv));
        Lnr = sum(log(nanray(origin, [dirxnr, dirynr, dirznr], rlim, ...
            lambda, xgv, ygv, zgv)));
        L(i,j) = Lr + Lnr;
        
        % Advance the progress bar.
        waitbar(((i-1)*numel(gvs) + j) / numel(gvs)^2, ...
            waitbarHandle);
    end
end
close(waitbarHandle);

%% Display result.
% Display the overall probabilities of the shifted scans.
surf(gvs, gvs, L);

% Add title and labels.
title('Log-likelihood of Lidar measurement')
xaxis('x [m]')
yaxis('y [m]')
