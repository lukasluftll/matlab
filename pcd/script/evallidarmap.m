% Compute KL divergence of lidar scans given a lidar map.

%% Set parameters.
% MAT file created by buildlidarmap script.
resultFolder = 'pcd/result';
lidarMapFile = [resultFolder, '/refmap_campus.mat'];

% Set the parameter that defines how many files make one scan.
pcdPerLs = 4;

% Minimum and maximum admissible map values.
mapLim = [0.002, 10];

%% Prepare output file.
% Load the file that contains the lidar map.
load(lidarMapFile);

% Define the name of the output MAT file.
[path, name, extension] = fileparts(lidarMapFile);
evalFile = [path, '/', name, '_eval', extension];

% Save parameters to file.
save(evalFile, 'lidarMapFile', 'pcdPerLs', 'mapLim', '-v7.3');

%% Compute KL divergence.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);

% Use appropriate sensor model.
switch model
    case 'decay'
        evalFun = @decayray;
    case 'ref'
        evalFun = @refray;
end

% Compute the KL divergence for each scan.
iScan = 1 : pcdPerLs : numel(pcdFile);
D = NaN(size(iScan));
parprogress(numel(iScan));
parfor i = 1 : numel(iScan)
    % Read laser scan data from files.
    ls = laserscan.empty(pcdPerLs, 0);
    for j = 1 : pcdPerLs
        filePath = [folder, '/', pcdFile(iScan(i)+j-1).name];
        ls(j) = lsread(filePath, rlim);
    end
    ls = lsconcat(ls);
    
    % Compute the measurement likelihood.
    [pi,Li] = evalFun(ls, constrain(lidarMap, mapLim));
    
    % Compute the KL divergence of the whole scan.
    D(i) = -sum([Li; log(pi)]);
    
    % Update the progress bar.
    parprogress;
end
parprogress(0);

% Save the KL divergence to file.
save(evalFile, 'D', '-append');
