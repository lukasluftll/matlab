% Visualize the log-likelihood of obtaining a shifted Lidar scan given a
% ray decay map created from the same Lidar scan at the true pose.
%
% Steps:
% # Load the scan.
% # Compute a Lidar decay map from it.
% # Shift the scan horizontally in x and y direction.
% # Compute the log-likelihood of obtaining the shifted scan with respect
%   to the decay map.

%% Set parameters.
% File name of the point cloud data file.
file = 'data/castle.pcd';

% Define the position of the sensor.
origin = [0, 0, 0];

% Shifting offset in x and y direction.
shift = 5;

% Resolution of the decay rate map.
res = 1;

% Resolution of the log-likelihood graph.
shiftres = 0.5;

% Minimum and maximum range of the Lidar sensor.
rlim = [0.9, 130];

% Maximum elevation angle of the Lidar sensor.
elevationMax = deg2rad(16);

% Minimum and maximum admissible decay rate.
lambdaLim = [2e-3, 1e+1];

%% Compute decay rate map.
% Read the point cloud.
pcd = pcdread(file);

% Compute the decay rate map.
radiusFinite = pcd.radius;
radiusFinite(~isfinite(radiusFinite)) = rlim(2);
hgv = -rlim(2)-shift : res : rlim(2)+res+shift;
vgv = -rlim(2)*sin(elevationMax) : res : rlim(2)*sin(elevationMax);
lambda = decaymap(origin, pcd.azimuth, pcd.elevation, radiusFinite, ...
    isfinite(pcd.radius), hgv, hgv, vgv);

% Set all voxels without data to the decay rate prior.
lambda.data(~isfinite(lambda.data)) = ...
    sum(isfinite(pcd.radius(:))) / sum(radiusFinite(:));

% Limit the decay rates to a reasonable interval.
lambda.data = max(lambdaLim(1), lambda.data);
lambda.data = min(lambdaLim(2), lambda.data);

% Visualize and save the decay rate map.
rayplot(pcd.azimuth, pcd.elevation, radiusFinite, isfinite(pcd.radius));
plot(lambda);
title('Decay rate map'); xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end
savefig(['pcd/results/decaymap_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);

%% Compute log-likelihood of shifted scans.
% Compute the direction vectors of the returned rays.
[dirxr, diryr, dirzr] = sph2cart(pcd.azimuth(isfinite(pcd.radius)), ...
    pcd.elevation(isfinite(pcd.radius)), pcd.radius(isfinite(pcd.radius)));

% Compute the direction vectors of the no-return rays.
[dirxnr, dirynr, dirznr] = sph2cart(pcd.azimuth(~isfinite(pcd.radius)), ...
    pcd.elevation(~isfinite(pcd.radius)), rlim(2));

% Shift the scan and compute the probability of obtaining it.
gvs = -shift : shiftres : shift;
L = zeros(numel(gvs));
for i = 1 : numel(gvs)
    for j = 1 : numel(gvs)
        origin = [gvs(i), gvs(j), 0];
        
        % Compute the log-likelihood of the reflected rays.
        Lr = sum(decayray(origin, [dirxr, diryr, dirzr], lambda));
        
        % Compute the log-likelihood of the no-return rays.
        Lnr = sum(log(decaynanray(origin, [dirxnr, dirynr, dirznr], ...
            rlim, lambda)));
        
        % Compute the log-likelihood of all measurements.
        L(i,j) = Lr + Lnr;
        
        % Display the overall probabilities of the shifted scans.
        surfHandle = surf(gvs, gvs, L);
        title('Log-likelihood of Lidar measurement from decay map')
        xlabel('x [m]')
        ylabel('y [m]')
        drawnow limitrate
    end
end

%% Save figure.
savefig(surfHandle, ['pcd/results/decaycastle_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);
