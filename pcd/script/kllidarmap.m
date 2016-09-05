% Compute KL divergence of lidar scans given a lidar map.

%% Fetch parameters.
lidarparams

%% Prepare output file.
% Load the file that contains the lidar map.
load(lidarMapFile, 'lidarMap', 'pnr');

% Print caption.
hline(75, '#')
display(['Computing KL divergence of ', model, 'map of ', dataset, ...
    ' dataset ...'])

%% Compute KL divergence.
try
    % Compute the KL divergence for each scan.
    Dgh = NaN(numel(evalFile), 1);
    parprogress(numel(evalFile));
    parfor i = 1 : numel(evalFile)
        % Read laser scan.
        ls = lsread([dataFolder, '/', evalFile(i).name], rlim);
        
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
    parprogress(0);
catch
end

% Save the KL divergence to file.
save(klFile, 'lidarMapFile', 'dataset', 'model', 'decayLim', ...
    'refLim', 'lfLim', 'Dgh', '-v7.3');
display(['Result written to ', klFile, '.'])
