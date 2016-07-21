% Visualizes th log-likelihood of obtaining a shifted Lidar scan given a
% reflectivity map created from the same Lidar scan at the true pose.
%
% Steps:
% # Load the scan.
% # Compute a Lidar reflectivity map from it.
% # Shift the scan horizontally in x and y direction.
% # Compute the log-likelihood of obtaining the shifted scan with respect
%   to the reflectivity map.

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

% Minimum and maximum admissible reflectivity.
refLim = [0.1, 0.9];

%% Compute decay rate map.
% Read the point cloud.
pcd = pcdread(file);

% Compute the decay rate map.
radiusFinite = pcd.radius;
radiusFinite(~isfinite(radiusFinite)) = rlim(2);
hgv = -rlim(2)-shift : res : rlim(2)+res+shift;
vgv = -rlim(2)*sin(elevationMax) : res : rlim(2)*sin(elevationMax);
ref = refmap(origin, pcd.azimuth, pcd.elevation, radiusFinite, ...
    isfinite(pcd.radius), hgv, hgv, vgv);

% Set all voxels without data to the reflectivity prior.
ref.data(~isfinite(ref.data)) = mean(ref.data(:), 'omitnan');

% Limit the decay rates to a reasonable interval.
ref.data = max(refLim(1), ref.data);
ref.data = min(refLim(2), ref.data);

% Visualize and save the reflectivity map.
plot(ref);
rayplot(pcd.azimuth, pcd.elevation, radiusFinite, isfinite(pcd.radius));
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end
savefig(['pcd/results/refmap_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);

%% Compute log-likelihood of shifted scans.
% Compute the direction vectors of the rays.
[dirx, diry, dirz] = sph2cart(pcd.azimuth, pcd.elevation, radiusFinite);

% Shift the scan and compute the probability of obtaining it.
gvs = -shift : shiftres : shift;
waitbarHandle = waitbar(0, 'Computing scan probabilities ...');
L = zeros(numel(gvs));
for i = 1 : numel(gvs)
    for j = 1 : numel(gvs)
        origin = [gvs(i), gvs(j), 0];
        
        % Compute the log-likelihood of the measurements.
        L(i,j) = sum(log(refray(origin, [dirx, diry, dirz], rlim, ref)));
        
        % Advance the progress bar.
        waitbar(((i-1)*numel(gvs) + j) / numel(gvs)^2, ...
            waitbarHandle);
    end
end
close(waitbarHandle);

%% Display result.
% Display the overall probabilities of the shifted scans.
surfHandle = surf(gvs, gvs, L);

% Add title and labels.
title('Log-likelihood of Lidar measurement from reflection map')
xlabel('x [m]')
ylabel('y [m]')

% Save figure.
savefig(['pcd/results/refcastle_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);
