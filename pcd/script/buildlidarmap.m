% Build a map out of many lidar scans.

%% Parameters.
% Define the folder from where to read and where to keep the results.
resultFolder = 'pcd/result';

% Name of the file that contains the merged point cloud.
pcFile = [resultFolder, '/pcmap_demo.mat'];

% Sensor model to use to build map: 'decay' | 'ref' | 'lf'.
model = 'ref';

% Resolution of the resulting lidar map.
mapRes = 5;

% Sensor reading range.
rlim = [2, 120];

% Variance for likelihood field computation.
sigma = 1;

%% Prepare output file.
% Load the file that contains the merged point cloud.
load(pcFile);

% Print caption.
hline(75)
disp(['Computing ', model, 'map for ', dataset, ' ...'])

% Define the name of the output MAT file.
lidarMapFile = [resultFolder, '/', model, 'map_', dataset, '.mat'];

%% Create lidar map.
% Compute the grid vectors of the map.
xgv = pcMap.XLimits(1)-pcRes : mapRes : pcMap.XLimits(2)+mapRes+pcRes;
ygv = pcMap.YLimits(1)-pcRes : mapRes : pcMap.YLimits(2)+mapRes+pcRes;
zgv = pcMap.ZLimits(2)-pcRes : mapRes : pcMap.ZLimits(2)+mapRes+pcRes;

% Define the mapping function to use.
switch lower(model)
    case 'decay'
        mapFun = @decaymap;
    case 'ref'
        mapFun = @refmap;
    case 'lf'
        mapFun = @lfmap;
    otherwise
        error('Mapping function not supported.')
end

% For the decay and the reflectivity model, loop over all scans and compute
% the numerator and denominator of each map cell's value.
% For the endpoint model, simply compute the Gaussian of the distance to
% the nearest obstacle for each voxel.
if strcmpi(model, 'decay') || strcmpi(model, 'ref')    
    % Get the PCD file names.
    pcdFile = dir([folder, '/*.pcd']);
    nPcdFile = numel(pcdFile);

    % Iterate over all laser scans, compute local maps and merge them into 
    % a global map.
    gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
    num = zeros(gridsize);
    denom = zeros(gridsize);
    parprogress(nPcdFile);
    parfor i = 1 : nPcdFile
        % Read laser scan data from file.
        ls = lsread([folder, '/', pcdFile(i).name], rlim);

        % Build the local lidar map.
        warning('off') %#ok<WNOFF>
        [~,numi,denomi] = mapFun(ls, xgv, ygv, zgv);
        warning('on') %#ok<WNON>

        % Integrate the local map information into the global map.
        num = num + numi;
        denom = denom + denomi;

        % Update progress bar.
        parprogress;
    end
    parprogress(0);
    
    % Compute the lidar map.
    lidarMap = voxelmap(num./denom, xgv, ygv, zgv);
else
    % Compute the likelihood field.
    lidarMap = mapFun(pcMap, sigma, xgv, ygv, zgv);
end

%% Save map.
save(lidarMapFile, 'dataset', 'folder', 'pcRes', 'model', ...
    'mapRes', 'rlim', 'sigma', 'lidarMap', '-v7.3');
display(['Result written to ', lidarMapFile, '.'])
