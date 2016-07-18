% Visualizes the log-likelihood of obtaining a shifted Lidar scan given a
% reflectivity map created from the same Lidar scan at the true pose.
%
% Steps:
% # Loads the scan.
% # Computes a Lidar reflectivity map from it.
% # Shifts the scan horizontally in x and y direction.
% # Computes the log-likelihood of obtaining the shifted scan with respect
%   to the reflectivity map.

%% Set parameters.
% Shifting offset in x and y direction.
shift = 3;

% Resolution of the decay rate map.
res = 0.5;

% Resolution of the log-likelihood graph.
shiftres = 0.5;

% Minimum and maximum range of the Lidar sensor.
rlim = [0, 130];

% Maximum elevation angle of the Lidar sensor.
elevationMax = deg2rad(16);

% Minimum and maximum admissible reflectivity.
refLim = [0, 1];

%% Compute decay rate map.
% Read the point cloud.
pcd = pcdread('data/castle.pcd');

% Compute the decay rate map.
radiusFinite = pcd.radius;
radiusFinite(~isfinite(radiusFinite)) = rlim(2);
hgv = -rlim(2)-shift : res : rlim(2)+res+shift;
zgv = -rlim(2)*sin(elevationMax) : res : rlim(2)*sin(elevationMax);
ref = rayref(pcd.azimuth, pcd.elevation, radiusFinite, ...
    isfinite(pcd.radius), hgv, hgv, zgv);

% Set all voxels without data to the reflectivity prior.
ref(~isfinite(ref)) = mean(ref(:));

% Limit the decay rates to a reasonable interval and add prior.
ref = max(refLim(1), ref);
ref = min(refLim(2), ref);

%% Compute log-likelihood of shifted scans.
% Compute the direction vectors of the returned rays.
[dirxr, diryr, dirzr] = sph2cart(pcd.azimuth(isfinite(pcd.radius)), ...
    pcd.elevation(isfinite(pcd.radius)), pcd.radius(isfinite(pcd.radius)));

% Compute the direction vectors of the no-return rays.
[dirxnr, dirynr, dirznr] = sph2cart(pcd.azimuth(~isfinite(pcd.radius)), ...
    pcd.elevation(~isfinite(pcd.radius)), rlim(2)*2);

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
        Lr = sum(pdfray(origin, [dirxr, diryr, dirzr], ref, ...
            hgv, hgv, zgv));
        Lnr = sum(log(nanray(origin, [dirxnr, dirynr, dirznr], rlim, ...
            ref, hgv, hgv, zgv)));
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
xlabel('x [m]')
ylabel('y [m]')
