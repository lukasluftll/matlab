% Visualize the log-likelihood of obtaining a shifted Lidar scan given a
% likelihood field created from the same Lidar scan at the true pose.
%
% Steps:
% # Load the scan.
% # Compute the likelihood field from it.
% # Shift the scan horizontally in x and y direction.
% # Compute the log-likelihood of obtaining the shifted scan with respect
%   to the likelihood field.

%% Set parameters.
% File name of the point cloud data file.
file = 'data/castle.pcd';

% Shifting offset in x and y direction.
shift = 5;

% Resolution of the decay rate map.
res = 1;

% Variance of the Gaussian kernel used to compute the likelihood field.
sigma = 1;

% Minimum and maximum range of the Lidar sensor.
rlim = [0.9, 130];

% Maximum elevation angle of the Lidar sensor.
elevationMax = deg2rad(16);

% Resolution of the log-likelihood graph.
shiftres = 0.5;

%% Compute decay rate map.
% Read the point cloud.
pcd = pcdread(file);
pc = pointCloud([pcd.x(:), pcd.y(:), pcd.z(:)]);

% Compute the likelihood field.
hgv = -rlim(2)-shift : res : rlim(2)+res+shift;
vgv = -rlim(2)*sin(elevationMax) : res : rlim(2)*sin(elevationMax);
lf = lfmap(pc, sigma, hgv, hgv, vgv);

% Visualize and save the likelihood field.
figure('Name', 'lfcastle map', 'NumberTitle', 'Off');
plot(laserscan(pcd.azimuth, pcd.elevation, pcd.radius));
plot(lf);
title('Likelihood field'); xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end
savefig(['pcd/results/lfmap_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);

%% Compute log-likelihood of shifted scans.
% Shift the scan and compute the probability of obtaining it.
gvs = -shift : shiftres : shift;
L = zeros(numel(gvs));
waitbarHandle = waitbar(0, 'Computing scan likelihoods ...');
for i = 1 : numel(gvs)
    for j = 1 : numel(gvs)
        % Compute the log-likelihood of the shifted point cloud.
        t = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; gvs(i), gvs(j), 0, 1]);
        L(i,j) = lfpc(pctransform(pc, t), lf);       
        
        % Advance the progress bar.
        waitbar(((i-1)*numel(gvs) + j) / numel(gvs)^2, waitbarHandle);
    end
end
close(waitbarHandle);

%% Plot log-likelihood of shifted scans.
% Display the overall log-likelihoods of the shifted scans.
figure('Name', 'lfcastle likelihood', 'NumberTitle', 'Off');
surf(gvs, gvs, L);
title('Log-likelihood of Lidar measurement from likelihood field')
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

% Save the figure.
savefig(['pcd/results/lfcastle_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);
