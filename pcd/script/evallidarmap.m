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

% Print caption.
hline(75)
display(['Evaluating ', model, 'map of ', dataset, ' dataset ...'])

% Define the name of the output MAT file.
[path, name, extension] = fileparts(lidarMapFile);
evalFile = [path, '/', name, '_eval', extension];

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
Dgh = NaN(size(iScan));
parprogress(numel(iScan));
parfor i = 1 : numel(iScan)
    % Read laser scan data from files.
    ls = laserscan.empty(pcdPerLs, 0);
    for j = 1 : pcdPerLs
        filePath = [folder, '/', pcdFile(iScan(i)+j-1).name]; %#ok<PFBNS>
        ls(j) = lsread(filePath, rlim);
    end
    ls = lsconcat(ls);
    
    % Compute the measurement likelihood.
    [pi,Li] = evalFun(ls, constrain(lidarMap, mapLim));
    
    % Compute the KL divergence of the whole scan.
    Dgh(i) = -sum([Li(~isnan(Li)); log(pi(~isnan(pi)))]);
    
    % Update the progress bar.
    parprogress;
end
parprogress(0);

% Save the KL divergence to file.
save(evalFile, 'lidarMapFile', 'pcdPerLs', 'mapLim', 'Dgh', '-v7.3');
display(['Result written to ', evalFile, '.'])
