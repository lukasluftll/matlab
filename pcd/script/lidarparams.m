% Parameters for scripts buildpcmap, buildlidarmap, evallidarmap.

%% Parameters.
% Resolution of the lidar map.
mapRes = 5;

% Resolution of the merged point cloud map.
pcRes = 0.1;

% Minimum and maximum admissible map values.
decayLim = [0.002, 10];
refLim = [0.001, 0.999];
lfLim = [1e-6, 10];

% Variance for likelihood field computation.
sigma = 1;

% Radius in x-y plane used for probability normalization when computing
% inverse KL divergence.
rkli = 2.5;

% Number of shifts used for probability normalization when computing
% inverse KL divergence.
nShift = 5;

% Ground-truth localization accuracy of the robot.
sigmaLoc = 0.2;

%% Derived parameters.
% Set the parameter that defines how many files make one scan.
pcdPerLs = 1;
switch dataset
    case 'demo'
        pcdPerLs = 1;
    case 'campus'
        pcdPerLs = 4;
    case 'schauinsland'
        pcdPerLs = 17;
    case 'mooswald'
        pcdPerLs = 19;
end

% Dataset folder containing the PCD files.
dataFolder = ['pcd/data/', dataset, '/pcd_sph'];

% Folder from where to read and where to keep the results.
resultFolder = 'pcd/result';

% Files used to build the map and localize the robot.
pcdFile = dir([dataFolder, '/*.pcd']);
iMap = rem(1:numel(pcdFile), 2*pcdPerLs) < pcdPerLs;
mappingFile = pcdFile(iMap);
evalFile = pcdFile(~iMap);

% Name of the MAT file that contains the merged point cloud.
pcMapFile = [resultFolder, '/pcmap_', dataset, '.mat'];

% Define default value for sensor model, if not defined.
if ~exist('model', 'var')
    model = 'decay';
end

% Name of the output file that contains the lidar map.
lidarMapFile = [resultFolder, '/', model, 'map_', dataset, '.mat'];

% Names of the evaluation files that contain KL divergence and inverse KL
% divergence.
klFile = [resultFolder, '/', model, 'map_', dataset, '_kl.mat'];
kliFile = [resultFolder, '/', model, 'map_', dataset, '_kli.mat'];

% Sensor reading range.
rlim = [2, 120];

%% Create folder for results.
[resultFolderPath, resultFolderName] = fileparts(resultFolder);
if ~exist(resultFolder, 'dir')
    mkdir(resultFolderPath, resultFolderName);
end
