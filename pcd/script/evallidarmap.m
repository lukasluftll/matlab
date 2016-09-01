% Compute KL divergence of lidar scans given a lidar map.

%% Prepare output file.
% Load the file that contains the lidar map.
load(lidarMapFile, 'lidarMap');

% Print caption.
hline(75)
display(['Evaluating ', model, 'map of ', dataset, ' dataset ...'])

% Define the name of the output MAT file.
[path, name, extension] = fileparts(lidarMapFile);
evalFile = [path, '/', name, '_eval', extension];

%% Compute KL divergence.
% Dataset folder containing the PCD files.
folder = ['pcd/data/', dataset, '/pcd_sph'];

% Get the PCD file names.
pcdFile = dir([folder, '/*.pcd']);

% Compute the KL divergence for each scan.
iScan = 1 : pcdPerLs : numel(pcdFile);
Dgh = NaN(size(iScan));
parprogress(numel(iScan));
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
parfor i = 1 : numel(iScan)
    % Read laser scan data from files.
    ls = laserscan.empty(pcdPerLs, 0); %#ok<*PFGV>
    for j = 1 : pcdPerLs
        filePath = [folder, '/', pcdFile(iScan(i)+j-1).name]; %#ok<PFBNS>
        ls(j) = lsread(filePath, rlim);
    end
    ls = lsconcat(ls);
    
    % Compute the measurement likelihood.
    switch lower(model)
        case 'decay'
            [pi, Li] = decayray(ls, constrain(lidarMap, decayLim));
        case 'ref'
            [pi, Li] = refray(ls, constrain(lidarMap, refLim));
        case 'lf'
            [pi, Li] = lfray(ls, constrain(lidarMap, lfLim), pnr);
        otherwise
            error(['Sensor model ', model, ' not supported.'])
    end
    
    % Compute the KL divergence of the whole scan.
    Dgh(i) = -sum([Li; log(pi)]);
    
    % Update the progress bar.
    parprogress;
end
warning('on', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
parprogress(0);

% Save the KL divergence to file.
save(evalFile, 'lidarMapFile', 'pcdPerLs', 'decayLim', 'refLim', ...
    'lfLim', 'Dgh', '-v7.3');
display(['Result written to ', evalFile, '.'])
