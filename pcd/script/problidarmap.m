% Compute probability distribution around selected position for different 
% sensor models.

%% Fetch parameters.
lidarparams

%% Prepare processing.
% Load the file that contains the lidar map.
load(lidarMapFile, 'lidarMap', 'pnr');

% Print caption.
hline(75, '#')
display(['Computing probability distribution for ', model, 'map of ', ...
    dataset, ' dataset ...'])

% Read a selected laser scan from file.
i = round(numel(evalFile) / 2);
ls = lsread([dataFolder, '/', evalFile(i).name], rlim);

%% Compute probability distribution.
ex = [];
try
    % Compute the probability of scans in grid around the true position.
    offset = 10;
    spacing = 1;
    xgv = -offset : spacing : offset;
    ygv = -offset : spacing : offset;
    parprogress(numel(xgv) * numel(ygv));
    L = NaN(numel(xgv), numel(ygv));
    parfor ix = 1 : numel(xgv)
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
            
            % Update the progress bar.
            parprogress;
        end
    end
    parprogress(0);
catch ex
    display(ex.msgtext)
end

% Save the KL divergence to file.
save(probFile, 'lidarMapFile', 'dataset', 'model', 'pcdPerLs', ...
    'decayLim', 'refLim', 'lfLim', 'L', 'xgv', 'ygv', 'ex', '-v7.3');
display(['Result written to ', probFile, '.'])
