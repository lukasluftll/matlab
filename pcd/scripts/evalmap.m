% Compute KL divergence of lidar scans given a lidar map.

%% Set parameters.
% MAT file created by buildmap script.
infile = 'pcd/results/decaymap_campus.mat';

% Step that determines the fraction of PCD files to use.
evalStep = 1;

% Minimum and maximum admissible reflection rate.
refLim = [0.001, 0.999];

% Minimum and maximum admissible decay rate.
lambdaLim = [0.002, 10];

%% Prepare output file.
% Load the input file.
load(infile);

% Create folder for results.
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end

% Define the name of the output MAT file.
[path, name, extension] = fileparts(infile);
outfile = [path, '/', name, '_eval', extension];

% Save parameters to file.
save(outfile, 'infile', 'evalStep', 'refLim', 'lambdaLim', '-v7.3');

%% Compute KL divergence.
% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);

% Create a progress bar.
waitbarHandle = waitbar(0, 'Computing KL divergence ...');

% Compute the KL divergence for each scan.
D = NaN(numel(pcdFile), 1);
for i = 1 : evalStep : numel(pcdFile)
    % Read laser scan data from file.
    ls = lsread([folder, '/', pcdFile(i).name], rlim);
    
    % Use the appropriate sensor model.
    switch model
        case 'decay'
            [pi,Li] = decayray(ls, constrain(lidarMap, lambdaLim));
        case 'ref'
            [pi,Li] = refray(ls, lidarMap);
    end
    
    % Compute the KL divergence of this lidar measurement.
    D(i) = -sum([Li; log(pi)]);
    
    % Save the KL divergence to file.
    save(outfile, 'D', '-append');
    
    % Advance the progress bar.
    waitbar(i/numel(pcdFile), waitbarHandle);
end

% Close the progress bar.
close(waitbarHandle);
