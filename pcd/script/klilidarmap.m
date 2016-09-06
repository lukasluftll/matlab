% Compute inverse KL divergence of lidar scans given a lidar map.

%% Fetch parameters.
lidarparams

%% Compute inverse KL divergence.
% Print caption.
hline(75, '#')
display(['Computing inverse KL divergence of ', model, 'map of ', ...
    dataset, ' dataset ...'])

% Load the lidar map and the unconditioned probability of no-return rays.
load(lidarMapFile, 'lidarMap', 'pnr');

% Catch any error that occurs while computing the inverse KL divergence.
ex = [];
try
    % Compute the inverse KL divergence for a subset of scans.
    n = numel(evalFile);
    fi = 1 : nShift*pcdPerLs : floor(n/pcdPerLs)*pcdPerLs;
    Dhg = NaN(numel(fi), 1);
    parprogress(numel(fi));    
    for i = 1 : numel(fi)
        % Read a full-revolution laser scan from file.
        ls = laserscan.empty(pcdPerLs, 0);
        for j = 1 : pcdPerLs
            ls(j) = lsread([dataFolder,'/',evalFile(fi(i)+j-1).name],rlim);
        end
        ls = lsconcat(ls);
        
        % Randomly shift the scan in the x-y plane and compute the overall
        % likelihood.
        j = 0;
        offset = NaN(nShift, 2);
        L = NaN(nShift, 1);
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
            warning('off', ...
                'MATLAB:mir_warning_maybe_uninitialized_temporary')
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
            warning('on', ...
                'MATLAB:mir_warning_maybe_uninitialized_temporary')
            
            % Compute the overall log-likelihood of obtaining the scan.
            L(j) = sum([log(pj); Lj]);
        end
        
        % Normalize the log-likelihood.
        L = L - logsumexp(L);
        
        % Normalize the Gaussian.
        N = mvnpdf(offset, 0, eye(2)*sigmaLoc);
        N = N / sum(N);
        
        % Compute the inverse KL divergence of the scan.
        Dhg(i) = sum(exp(L) .* (L - log(N)));
        
        % Compute baseline inverse KL divergence resulting from a uniform
        % distribution.
        Dhgbl = -sum(log(nShift) + log(N)) / nShift;
        
        % Update the progress bar.
        parprogress;
    end
    parprogress(0);
catch ex
    display(ex.message)
end

% Save the KL divergence to file.
save(kliFile, 'lidarMapFile', 'dataset', 'model', 'pcdPerLs', ...
    'decayLim', 'refLim', 'lfLim', 'Dhg', 'Dhgbl', 'ex', '-v7.3');
display(['Result written to ', kliFile, '.'])
