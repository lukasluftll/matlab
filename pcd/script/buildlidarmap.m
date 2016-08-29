% Build a map out of many lidar scans.

%% Parameters.
% Define the folder from where to read and where to keep the results.
resultFolder = 'pcd/result';

% Name of the file that contains the merged point cloud.
pcFile = [resultFolder, '/pcmap_demo.mat'];

% Sensor model to use to build map: 'decay' | 'ref' | 'lf'.
model = 'decay';

% Resolution of the resulting lidar map.
mapRes = 2.0;

% Sensor reading range.
rlim = [2, 120];

% Variance for likelihood field computation.
sigma = 1;

%% Prepare output file.
% Load the file that contains the merged point cloud.
load(pcFile);

% Print caption.
hline(75)
disp(['Computing ', model, 'map for ', dataset, ' dataset ...'])

% Define the name of the output MAT file.
lidarMapFile = [resultFolder, '/', model, 'map_', dataset, '.mat'];

%% Compute extent of lidar map.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);
nPcdFile = numel(pcdFile);

% Iterate over all scans and compute an axis-aligned bounding box of all
% sensor poses.
mapLim = repmat([Inf, -Inf], 3, 1);
for i = 1 : nPcdFile
    ls = lsread([folder, '/', pcdFile(i).name], rlim);
    spmin = min(tform2trvec(ls.sp)).';
    spmax = max(tform2trvec(ls.sp)).';
    mapLim = [min([mapLim(:,1), spmin], [], 2), ...
        max([mapLim(:,2), spmax], [], 2)];
end

% Compute the bounding box of a scan that consists of maximum range
% readings only.
elevationLim = [min(ls.elevation), max(ls.elevation)];
zlim = rlim(2) * sin(elevationLim);
xylim = [-1,+1] * rlim(2);

% Add the extent of the maximum range scan to the bounding box of the
% sensor poses to get the extent of the map.
mapLim = mapLim + [xylim; xylim; zlim];

%% Create lidar map.
% Compute the grid vectors of the map.
xgv = mapLim(1,1) : mapRes : mapLim(1,2)+mapRes;
ygv = mapLim(2,1) : mapRes : mapLim(2,2)+mapRes;
zgv = mapLim(3,1) : mapRes : mapLim(3,2)+mapRes;

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
