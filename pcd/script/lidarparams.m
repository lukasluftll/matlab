% Parameters for scripts buildpcmap, buildlidarmap, evallidarmap.

%% Dynamic parameters.
% Dataset name.
dataset = 'demo';

% Sensor model to use to build the map: 'decay' | 'ref' | 'lf'.
model = 'decay';

% Set the parameter that defines how many files make one scan.
pcdPerLs = 4;

% Resolution of the lidar map.
mapRes = 0.5;

% Resolution of the merged point cloud map.
pcRes = 0.1;

% Minimum and maximum admissible map values.
decayLim = [0.002, 10];
refLim = [0.001, 0.999];
lfLim = [1e-6, 10];

% Variance for likelihood field computation.
sigma = 1;

%% Static parameters.
% Folder from where to read and where to keep the results.
resultFolder = 'pcd/result';

% Sensor reading range.
rlim = [2, 120];
