% Build a map out of many lidar scans.

%% Parameters.
% Dataset name.
dataset = 'campus';

% Sensor model to use to build map: 'decay' | 'ref'.
model = 'decay';

% Dataset folder with PCD files.
folder = ['pcd/data/', dataset, '/pcd_sph'];

% Step that determines the fraction of PCD files to use.
step = 1;

% Resolution of the merged point cloud map.
pcMapRes = 0.5;

% Resolution of the resulting lidar map.
lidarMapRes = 0.5;

% Sensor reading range.
rlim = [2, 120];

%% Validate parameters.
switch model
    case 'decay'
    case 'ref'
    otherwise
        error(['Map model ', model, ' not supported.'])
end

%% Prepare output file.
% Create folder for results.
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end

% Define the name of the output MAT file.
outfile = ['pcd/results/', model, 'map_', dataset, '.mat'];

% Save parameters to file.
save(outfile, 'dataset', 'model', 'folder', 'step', 'pcMapRes', ...
    'lidarMapRes', 'rlim', '-v7.3');

%% Merge point clouds.
% Get the PCD file names.
infile = dir([folder, '/*.pcd']);

% Create progress bar.
waitbarHandle = waitbar(0, 'Merging point cloud map ...');

% Iterate over all PCD files.
pcMap = pointCloud(zeros(0, 3));
for i = 1 : step : numel(infile)
    % Read laser scan data from file.
    ls = lsread([folder, '/', infile(i).name], rlim);
    
    % Merge this point cloud with the map.
    pcMap = pcmerge(pcMap, ls2pc(ls), pcMapRes);
    
    % Advance the progress bar.
    waitbar(i/numel(infile), waitbarHandle);
end

% Save the point cloud map to file.
save(outfile, 'pcMap', '-append');

% Denoise the map.
pcMap = pcdenoise(pcMap);

% Compute the extent of the denoised point cloud.
lim = [pcMap.XLimits(1), pcMap.XLimits(2);
    pcMap.YLimits(1), pcMap.YLimits(2);
    pcMap.ZLimits(1), pcMap.ZLimits(2)];

%% Create lidar map.
% Compute the grid vectors of the map.
xgv = lim(1,1) : lidarMapRes : lim(1,2)+lidarMapRes;
ygv = lim(2,1) : lidarMapRes : lim(2,2)+lidarMapRes;
zgv = lim(3,1) : lidarMapRes : lim(3,2)+lidarMapRes;

% Update progress bar label.
waitbar(0, waitbarHandle, ['Computing ', model, ' map ...']);

% Create the lidar map.
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
numerator = voxelmap(zeros(gridsize), xgv, ygv, zgv);
denominator = numerator.copy;
for i = 1 : step : numel(infile)
    % Read laser scan data from file.
    ls = lsread([folder, '/', infile(i).name], rlim);
    
    % Build the local lidar map.
    switch model
        case 'decay'
            [~,ai,bi] = decaymap(ls, xgv, ygv, zgv);
        case 'ref'
            [~,ai,bi] = refmap(ls, xgv, ygv, zgv);
    end
    
    % Integrate the local map information into the global map.
    numerator.add(ai);
    denominator.add(bi);
    lidarMap = numerator ./ denominator;
    
    % Save the global map to file.
    save(outfile, 'numerator', 'denominator', 'lidarMap', '-append');
    
    % Advance the progress bar.
    waitbar(i/numel(infile), waitbarHandle);
end

% Close the progress bar.
close(waitbarHandle);
