% Compute inverse KL divergence of lidar scans given a lidar map.

%% Fetch parameters.
lidarparams

%% Prepare output file.
% Load the file that contains the lidar map.
load(lidarMapFile, 'lidarMap');

% Print caption.
hline(75, '#')
display(['Computing inverse KL divergence of ', model, 'map of ', ...
    dataset, ' dataset ...'])

% Define the name of the output MAT file.
[path, name, extension] = fileparts(lidarMapFile);
evalFile = [path, '/', name, '_kli', extension];

%% Compute KL divergence.
% Get the PCD file names.
pcdFile = dir([folder, '/evaluation/*.pcd']);

% Compute the inverse KL divergence for each scan.
nPcdFile = numel(pcdFile);
Dhg = NaN(nPcdFile, 1);
parprogress(nPcdFile);
warning('off', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
for i = 1 : nPcdFile
    % Read laser scan data from file.
    ls = lsread([folder, '/evaluation/', pcdFile(i).name], rlim);
    
    % Randomly shift the scan in the x-y plane and compute the overall
    % likelihood.
    j = 0;
    offset = NaN(nShift, 2);
    ps = NaN(nShift, 1);
    while j < nShift
        % Generate random x-y offset.
        offsetTmp = rkli * ([2*rand,2*rand]-1);
        
        % Make sure offset position is within radius.
        if norm(offsetTmp) > rkli
            continue
        end
        
        % Increment number of shifts.
        j = j + 1;
        
        % Shift laser scan.
        offset(j,:) = offsetTmp;
        lss = ls;
        lss.sp = pagetimes(lss.sp, ...
            repmat(trvec2tform([offset(j,:),0]), 1, 1, size(lss.sp,3)));
        
        % Compute the measurement likelihood.
        switch lower(model)
            case 'decay'
                [pr, Lr] = decayray(lss, constrain(lidarMap, decayLim));
            case 'ref'
                [pr, Lr] = refray(lss, constrain(lidarMap, refLim));
            case 'lf'
                [pr, Lr] = lfray(lss, constrain(lidarMap, lfLim), pnr);
            otherwise
                error(['Sensor model ', model, ' not supported.'])
        end
        
        % Sum up the overall probability of the scan.
        Ls = sum([log(pr); Lr]);
        ps(j) = exp(Ls);
    end
    
    % Normalize the probabilities.
    ptot = sum(ps);
    ps = ps / ptot;
    
    % Normalize the Gaussian.
    Ntot = sum(mvnpdf(offset, [0,0], eye(2)*sigma));
    Ns = Ns / Ntot;

    % Compute the inverse KL divergence of the whole scan.
    Dhg(i) = sum(ps * (log(ps) - log(Ns)));
    
    % Update the progress bar.
    parprogress;
end
warning('on', 'MATLAB:mir_warning_maybe_uninitialized_temporary')
parprogress(0);

% Save the KL divergence to file.
save(evalFile, 'lidarMapFile', 'pcdPerLs', 'decayLim', 'refLim', ...
    'lfLim', 'Dgh', '-v7.3');
display(['Result written to ', evalFile, '.'])
