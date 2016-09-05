% Compute KL divergence of lidar scans given a lidar map.

%% Fetch parameters.
lidarparams

%% Compute KL divergence.
% Print caption.
hline(75, '#')
display(['Computing KL divergence of ', model, 'map of ', dataset, ...
    ' dataset ...'])

% Load the lidar map and the unconditioned probability of no-return rays.
load(lidarMapFile, 'lidarMap', 'pnr');

% Determine the number of files used for evaluation.
n = numel(evalFile);

% Catch any error that occurs while computing the KL divergence.
ex = [];
try
    % Compute the KL divergence for each scan.
    parprogress(n);
    Dghs = NaN(n, 1);
    parfor i = 1 : n
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
        
        % Compute the KL divergence of the scan.
        Dghs(i) = -sum([log(pi); Li]);
        
        % Update the progress bar.
        parprogress;
    end
    parprogress(0);
catch ex
    display(ex.message)
end

% Merge the KL divergences of several scans to get the divergence of a 
% full scanner revolution.
Dgh = sum(reshape(Dghs(1 : floor(n/pcdPerLs)*pcdPerLs), pcdPerLs, []));

% Save the KL divergence to file.
save(klFile, 'lidarMapFile', 'dataset', 'model', 'decayLim', ...
    'refLim', 'lfLim', 'pcdPerLs', 'Dghs', 'Dgh', 'ex', '-v7.3');
display(['Result written to ', klFile, '.'])
