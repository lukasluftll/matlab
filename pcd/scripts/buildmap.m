% Build a map out of many lidar scans.

%% Set parameters.
% Dataset name.
dataset = 'campus';

% Sensor model to use to build map: 'decay' | 'ref'.
model = 'decay';

% Step that determines the fraction of PCD files to use.
step = 1;

% Resolution of the merged point cloud map.
pcres = 0.200;

% Resolution of the resulting lidar map.
lambdares = 0.200;

% Sensor reading range.
rlim = [2, 120];

% Save parameters to file.
save(['pcd/results/', model, 'map_', dataset, '.mat'], ...
    'dataset', 'model', 'step', 'pcres', 'lambdares', 'rlim');

%% Create folder for results.
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end

%% Create one point cloud of all scans.
% Get the PCD file names.
folder = ['pcd/data/', dataset, '/pcd_sph'];
file = dir([folder, '/*.pcd']);

% Create progress bar.
waitbarHandle = waitbar(0, 'Building point cloud map ...');

% Iterate over all PCD files.
pcmap = pointCloud(zeros(0, 3));
for i = 1 : step : numel(file)
    % Read laser scan data from file.
    ls = lsread([folder, '/', file(i).name], rlim);
    
    % Merge this point cloud with the map.
    pcmap = pcmerge(pcmap, ls2pc(ls), pcres);
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

% Save the point cloud map to file.
save(['pcd/results/', model, 'map_', dataset, '.mat'], 'pcmap');

% Denoise the map.
pcmap = pcdenoise(pcmap);

% Compute the extent of the denoised point cloud.
lim = [pcmap.XLimits(1), pcmap.XLimits(2);
    pcmap.YLimits(1), pcmap.YLimits(2);
    pcmap.ZLimits(1), pcmap.ZLimits(2)];

%% Create lidar map.
% Compute the grid vectors of the map.
xgv = lim(1,1) : lambdares : lim(1,2)+lambdares;
ygv = lim(2,1) : lambdares : lim(2,2)+lambdares;
zgv = lim(3,1) : lambdares : lim(3,2)+lambdares;

% Update progress bar label.
waitbar(0, waitbarHandle, ['Building ', model, ' map ...']);

% Create the lidar map.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
a = voxelmap(zeros(gridsize), xgv, ygv, zgv);
b = a.copy;
for i = 1 : step : numel(file)
    % Read laser scan data from file.
    ls = lsread([folder, '/', file(i).name], rlim);
    
    % Build the local lidar map.
    switch model
        case 'decay'
            [~,ai,bi] = decaymap(ls, xgv, ygv, zgv);
        case 'ref'
            [~,ai,bi] = refmap(ls, xgv, ygv, zgv);
        otherwise
            error(['Map model ', model, ' not supported.'])
    end
    
    % Integrate the local map information into the global map.
    a.add(ai);
    b.add(bi);
    
    % Compute the global decay rate map.
    lidarmap = a ./ b;
    
    % Save the decay rate map to file.
    save(['pcd/results/', model, 'map_', dataset, '.mat'], 'lidarmap', ...
        '-append');
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

% Close the progress bar.
close(waitbarHandle);
