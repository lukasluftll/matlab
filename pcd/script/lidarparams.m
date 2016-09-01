% Parameters for scripts buildpcmap, buildlidarmap, evallidarmap.

%% Dynamic parameters.
% Dataset name.
dataset = 'demo';

% Sensor model to use to build the map: 'decay' | 'ref' | 'lf'.
model = 'decay';

% Set the parameter that defines how many files make one scan.
pcdPerLs = 1;

% Resolution of the lidar map.
mapRes = 1;

% Resolution of the merged point cloud map.
pcRes = 0.1;

% Minimum and maximum admissible map values.
decayLim = [0.002, 10];
refLim = [0.001, 0.999];
lfLim = [1e-6, 10];

% Variance for likelihood field computation.
sigma = 1;

%% Static parameters.
% Dataset folder containing the PCD files.
folder = ['pcd/data/', dataset, '/pcd_sph'];

% Folder from where to read and where to keep the results.
resultFolder = 'pcd/result';

% Sensor reading range.
rlim = [2, 120];

% MAT file that contains the merged point cloud.
pcFile = [resultFolder, '/pcmap_', dataset, '.mat'];

% MAT file that contains the lidar map.
lidarMapFile = [resultFolder, '/', model, 'map_', dataset, '.mat'];
