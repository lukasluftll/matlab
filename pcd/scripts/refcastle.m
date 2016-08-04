% Visualizes the log-likelihood of obtaining a shifted Lidar scan given a
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

% Define the position of the sensor w.r.t. the reflectivity map.
mts = affine3d();

% Shifting offset in x and y direction.
shift = 5;

% Resolution of the reflectivity map.
res = 1;

% Resolution of the log-likelihood graph.
shiftres = 0.5;

% Minimum and maximum range of the Lidar sensor.
rlim = [0.9, 130];

% Maximum elevation angle of the Lidar sensor.
elevationMax = deg2rad(16);

% Minimum and maximum admissible reflectivity.
refLim = [0.001, 0.999];

%% Compute reflectivity map.
% Read the point cloud.
pcd = pcdread(file);

% Compute the reflectivity map.
radiusFinite = pcd.radius;
radiusFinite(~isfinite(pcd.radius)) = rlim(2);
hgv = -rlim(2)-shift : res : rlim(2)+res+shift;
vgv = -rlim(2)*sin(elevationMax) : res : rlim(2)*sin(elevationMax);
ref = refmap(mts, pcd.azimuth, pcd.elevation, radiusFinite, ...
    isfinite(pcd.radius), hgv, hgv, vgv);

% Set all voxels without data to the reflectivity prior.
ref.data(~isfinite(ref.data)) = mean(ref.data(:), 'omitnan');

% Limit the decay rates to the admissible reflectivity interval.
ref.data = constrain(ref.data, refLim);

% Visualize and save the reflectivity map.
figure('Name', 'refcastle map', 'NumberTitle', 'Off')
rayplot(pcd.azimuth, pcd.elevation, radiusFinite, isfinite(pcd.radius));
plot(ref);
title('Reflectivity map'); xlabel('x[m]'); ylabel('y[m]'); zlabel('z[m]');
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end
savefig(['pcd/results/refmap_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);

%% Compute log-likelihood of shifted scans.
% Set the length of the NaN rays to a value outside the measurement range.
radiusFinite(~isfinite(pcd.radius)) = 1.1 * rlim(2);

% Compute the direction vectors of the rays.
[dirx, diry, dirz] = sph2cart(pcd.azimuth, pcd.elevation, radiusFinite);

% Shift the scan and compute the probability of obtaining it.
gvs = -shift : shiftres : shift;
L = zeros(numel(gvs));
waitbarHandle = waitbar(0, 'Computing scan likelihoods ...');
for i = 1 : numel(gvs)
    for j = 1 : numel(gvs)
        origin = [gvs(i), gvs(j), 0];
        
        % Compute the log-likelihood of the measurements.
        L(i,j) = sum(log(refray(origin, [dirx, diry, dirz], rlim, ref)));
        
        % Advance progress bar.
        waitbar(((i-1)*numel(gvs) + j) / numel(gvs)^2, waitbarHandle);
    end
end
close(waitbarHandle);

%% Plot log-likelihood of shifted scans.
% Display the overall log-likelihoods of the shifted scans.
figure('Name', 'refcastle likelihood', 'NumberTitle', 'Off');
surf(gvs, gvs, L);
title('Log-likelihood of Lidar measurement from reflection map')
xlabel('x [m]'); ylabel('y [m]'); zlabel('z [m]');

% Save the figure.
savefig(['pcd/results/refcastle_', ...
    datestr(now, 'yyyy-mm-dd_HH-MM-SS'), '.fig']);
