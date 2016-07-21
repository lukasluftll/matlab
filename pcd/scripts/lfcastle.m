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
rayplot(pcd.azimuth, pcd.elevation, pcd.radius);
plot(lf);
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end
savefig(['pcd/results/lfmap_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);

%% Compute log-likelihood of shifted scans.
% Shift the scan and compute the probability of obtaining it.
gvs = -shift : shiftres : shift;
L = zeros(numel(gvs));
for i = 1 : numel(gvs)
    for j = 1 : numel(gvs)
        % Shift the point cloud horizontally.
        t = affine3d([1 0 0 0; 0 1 0 0; 0 0 1 0; gvs(i), gvs(j), 0, 1]);
        pct = pctransform(pc, t);
        
        % Compute the log-likelihood of all measurements.
        L(i,j) = lfpc(pct, lf);
    end
end

%% Display result.
% Display the overall probabilities of the shifted scans.
surfHandle = surf(gvs, gvs, L);

% Add title and labels.
title('Log-likelihood of Lidar measurement from likelihood field')
xlabel('x [m]')
ylabel('y [m]')

% Save figure.
savefig(['pcd/results/lfcastle_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);
