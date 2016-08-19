% Compute KL divergence of lidar scans given a lidar map.

%% Set parameters.
% Dataset name.
dataset = 'campus';

% Sensor model: 'decay' | 'ref'.
model = 'decay';

% Step that determines the fraction of PCD files to use.
step = 1000;

% Sensor reading range.
rlim = [2, 120];

% Minimum and maximum admissible decay rate.
lambdaLim = [2e-3, 1e+1];

%% Create folder for results.
if ~exist('pcd/results', 'dir')
    mkdir('pcd', 'results');
end

%% Compute KL divergence.
% Get the PCD file names.
folder = ['pcd/data/', dataset, '/pcd_sph'];
file = dir([folder, '/*.pcd']);

% Create a progress bar.
waitbarHandle = waitbar(0, 'Computing KL divergence ...');

% Compute the KL divergence for each scan.
D = [];
for i = 1 : step : numel(file)
    % Read laser scan data from file.
    ls = lsread([folder, '/', file(i).name], rlim);
    
    % Use the appropriate sensor model.
    switch model
        case 'decay'
            [pi,Li] = decayray(ls, constrain(lidarmap, lambdaLim));
        case 'ref'
            [pi,Li] = refray(ls, lidarmap);
        otherwise
            error(['Map model ', model, ' not supported.'])
    end
    
    % Compute the KL divergence of this lidar measurement.
    D = [D; sum([Li; log(pi)])]; %#ok<AGROW>
    
    % Save the KL divergence to file.
    save(['pcd/results/kl', model, '_', dataset, '.mat'], 'D');
    
    % Advance the progress bar.
    waitbar(i/numel(file), waitbarHandle);
end

% Close the progress bar.
close(waitbarHandle);
