% Compute inverse KL divergence of lidar scans given a lidar map.

%% Fetch parameters.
lidarparams

%% Prepare processing.
% Load the file that contains the lidar map.
load(lidarMapFile, 'lidarMap', 'pnr');

% Print caption.
hline(75, '#')
display(['Computing inverse KL divergence of ', model, 'map of ', ...
    dataset, ' dataset ...'])

%% Compute inverse KL divergence.
try
    % Compute the inverse KL divergence for each scan.
    fi = 1 : nShift : numel(evalFile);
    Dhg = NaN(numel(fi), 1);
    parprogress(numel(fi));
    for i = 1 : numel(fi)
        % Read laser scan data from file.
        ls = lsread([dataFolder, '/', evalFile(fi(i)).name], rlim);
        
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
                repmat(trvec2tform([offset(j,:),0]),1,1,size(lss.sp,3)));
            
            % Compute the measurement likelihood.
            switch lower(model)
                case 'decay'
                    [pj, Lj] = decayray(lss, constrain(lidarMap,decayLim));
                case 'ref'
                    [pj, Lj] = refray(lss, constrain(lidarMap,refLim));
                case 'lf'
                    [pj, Lj] = lfray(lss, constrain(lidarMap,lfLim), pnr);
                otherwise
                    error(['Sensor model ', model, ' not supported.'])
            end
            
            % Integrate the probability of the scan over the circle area.
            ps(j) = sum([pj; exp(Lj)]);
        end
        
        % Normalize the probabilities.
        ps = ps / sum(ps);
        
        % Normalize the Gaussian.
        Ns = mvnpdf(offset, 0, eye(2)*sigmaLoc);
        Ns = Ns / sum(Ns);
        
        % Compute the inverse KL divergence of the scan.
        Dhg(i) = sum(ps .* (log(ps) - log(Ns)));
        
        % Update the progress bar.
        parprogress;
    end
    parprogress(0);
catch
end

% Save the KL divergence to file.
save(kliFile, 'lidarMapFile', 'dataset', 'model', 'pcdPerLs', ...
    'decayLim', 'refLim', 'lfLim', 'Dhg', '-v7.3');
display(['Result written to ', kliFile, '.'])
