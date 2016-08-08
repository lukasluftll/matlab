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
ls = laserscan(pcd.azimuth, pcd.elevation, pcd.radius, rlim);
hgv = -rlim(2)-shift : res : rlim(2)+res+shift;
vgv = -rlim(2)*sin(elevationMax) : res : rlim(2)*sin(elevationMax);
lambda = decaymap(ls, hgv, hgv, vgv);

% Set all voxels without data to the decay rate prior.
lambda.data(~isfinite(lambda.data)) = ...
    sum(isfinite(pcd.radius(:))) / sum(radiusFinite(:));

% Limit the decay rates to a reasonable interval.
lambda.data = max(lambdaLim(1), lambda.data);
lambda.data = min(lambdaLim(2), lambda.data);

% Visualize and save the decay rate map.
figure('Name', 'decaycastle map', 'NumberTitle', 'Off');
plot(ls);
plot(lambda);
title('Decay rate map'); xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end
savefig(['pcd/results/decaymap_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);

%% Compute log-likelihood of shifted scans.
% Shift the scan and compute the probability of obtaining it.
gvs = -shift : shiftres : shift;
L = zeros(numel(gvs));
waitbarHandle = waitbar(0, 'Computing scan likelihoods ...');
for i = 1 : numel(gvs)
    for j = 1 : numel(gvs)
        origin = [gvs(i), gvs(j), 0];
        
        % Compute the log-likelihood of the measurements.
        ls.tform(1:2,4) = [gvs(i); gvs(j)];
        [Ls, ps] = decayray(ls, lambda);
              
        % Compute the log-likelihood of all measurements.
        L(i,j) = sum([Ls; log(ps)]);
        
        % Advance the progress bar.
        waitbar(((i-1)*numel(gvs) + j) / numel(gvs)^2, waitbarHandle);
    end
end
close(waitbarHandle);

%% Plot log-likelihood of shifted scans.
% Display the overall log-likelihoods of the shifted scans.
figure('Name', 'decaycastle likelihood', 'NumberTitle', 'Off');
surf(gvs, gvs, L);
title('Log-likelihood of Lidar measurement from decay map');
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

% Save the figure.
savefig(['pcd/results/decaycastle_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);
