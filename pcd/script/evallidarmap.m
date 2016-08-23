% Compute KL divergence of lidar scans given a lidar map.

%% Set parameters.
% MAT file created by buildmap script.
resultFolder = 'pcd/result';
lidarMapFile = [resultFolder, '/decaymap_campus.mat'];

% Set the step that determines the fraction of PCD files to use.
evalStep = 1;

% Minimum and maximum admissible map values.
mapLim = [0.002, 10];

%% Prepare output file.
% Load the input file.
load(lidarMapFile);

% Define the name of the output MAT file.
[path, name, extension] = fileparts(lidarMapFile);
evalFile = [path, '/', name, '_eval', extension];

% Save parameters to file.
save(evalFile, 'lidarMapFile', 'evalStep', 'mapLim', '-v7.3');

%% Compute KL divergence.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);

% Create a progress bar.
progressbar('Computing KL divergence ...');

% Use appropriate sensor model.
switch model
    case 'decay'
        evalFun = @decayray;
    case 'ref'
        evalFun = @refray;
end

% Compute the KL divergence for each scan.
D = NaN(numel(pcdFile), 1);
nPcdFile = numel(pcdFile);
parfor i = 1 : evalStep : nPcdFile
    % Read laser scan data from file.
    ls = lsread([folder, '/', pcdFile(i).name], rlim);
    
    % Compute the measurement likelihood.
    [pi,Li] = evalFun(ls, constrain(lidarMap, mapLim));
    
    % Compute the KL divergence of this lidar measurement.
    D(i) = -sum([Li; log(pi)]);
    
    % Advance the progress bar.
    progressbar(i/nPcdFile);
end

% Save the KL divergence to file.
save(evalFile, 'D', '-append');
