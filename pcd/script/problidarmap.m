% Compute probability distribution around selected position for different 
% sensor models.

%% Fetch parameters.
lidarparams

%% Compute probability distribution.
% Print caption.
hline(75, '#')
display(['Computing probability distribution for ', model, 'map of ', ...
    dataset, ' dataset ...'])

% Load the lidar map and the unconditioned probability of no-return rays.
load(lidarMapFile, 'lidarMap', 'pnr');

% Read a selected full-revolution laser scan from file.
n = numel(evalFile);
fi = 1 : nShift*pcdPerLs : floor(n/pcdPerLs)*pcdPerLs;
fi = fi(ceil(numel(fi)/2));
ls = laserscan.empty(pcdPerLs, 0);
for i = 1 : pcdPerLs
    ls(i) = lsread([dataFolder, '/', evalFile(fi+i-1).name], rlim);
end
ls = lsconcat(ls);
        
% Catch any error.
ex = [];
try
    % Compute the log-likelihood of scans in a grid around the true 
    % position.
    spacing = 1;
    xgv = -rkli : spacing : rkli;
    ygv = -rkli : spacing : rkli;
    parprogress(numel(xgv) * numel(ygv));
    L = NaN(numel(xgv), numel(ygv));
    for ix = 1 : numel(xgv)
        for iy = 1 : numel(ygv)
            % Shift laser scan.
            lss = ls;
            lss.sp = pagetimes(lss.sp, repmat(trvec2tform(...
                [xgv(ix),ygv(iy),0]),1,1,size(lss.sp,3)));
            
            % Compute the measurement likelihood.
            switch lower(model)
                case 'decay'
                    [pi,Li] = decayray(lss, constrain(lidarMap,decayLim));
                case 'ref'
                    [pi,Li] = refray(lss, constrain(lidarMap,refLim));
                case 'lf'
                    [pi,Li] = lfray(lss, constrain(lidarMap,lfLim), pnr);
                otherwise
                    error(['Sensor model ', model, ' not supported.'])
            end
            
            % Compute the log-likelihood of the whole scan.
            L(ix,iy) = sum([log(pi); Li]);
            
            parprogress;
        end
    end
    parprogress(0);
catch ex
    display(ex.message)
end

% Save the KL divergence to file.
save(probFile, 'lidarMapFile', 'dataset', 'model', 'pcdPerLs', ...
    'decayLim', 'refLim', 'lfLim', 'L', 'xgv', 'ygv', 'ex', '-v7.3');
display(['Result written to ', probFile, '.'])

%% Visualize results.
% To display the log-likelihood, use:
% surf(xgv, ygv, L.')
