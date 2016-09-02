% Compute KL divergence of lidar scans given a lidar map.

%% Fetch parameters.
lidarparams

%% Prepare output file.
% Load the file that contains the lidar map.
load(lidarMapFile, 'lidarMap');

% Print caption.
hline(75, '#')
display(['Computing KL divergence of ', model, 'map of ', dataset, ...
    ' dataset ...'])

%% Compute KL divergence.
% Get the PCD file names.
pcdFile = dir([evaluationFolder, '/*.pcd']);

% Compute the KL divergence for each scan.
iScan = 1 : pcdPerLs : numel(pcdFile);
Dgh = NaN(size(iScan));
parprogress(numel(iScan));
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
parfor i = 1 : numel(iScan)
    % Read multiple laser scans from files and concatenate to get a full
    % revolution.
    ls = laserscan.empty(pcdPerLs, 0); %#ok<*PFGV>
    for j = 1 : pcdPerLs %#ok<*PFBNS>
        filePath = [evaluationFolder, '/', pcdFile(iScan(i)+j-1).name]; 
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
save(evalKlFile, 'lidarMapFile', 'pcdPerLs', 'decayLim', 'refLim', ...
    'lfLim', 'Dgh', '-v7.3');
display(['Result written to ', evalKlFile, '.'])
