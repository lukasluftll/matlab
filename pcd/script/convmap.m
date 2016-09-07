% Check convergence of a subvolume of the campus map.

%% Set parameters.
dataset = 'campus';
lidarparams
mapRes = 1;
xgv = -40 : mapRes : 0;
ygv =  10 : mapRes : 40;
zgv =  -5 : mapRes : 5;

%% Compute map convergence.
% Print caption.
hline(75, '#')
disp(['Computing ',model,'map convergence for ',dataset,' dataset ...'])

%% Create lidar maps.
% Loop over all scans and compute local maps.
n = numel(mappingFile);
gridsize = [numel(xgv), numel(ygv), numel(zgv)] - 1;
r = zeros([gridsize, n]);
l = zeros([gridsize, n]);
h = zeros([gridsize, n]);
m = zeros([gridsize, n]);
parprogress(n);
parfor i = 1 : n
    % Read laser scan data from file.
    ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);

    % Build the local decay rate map.
    warning('off', 'pcd:mapping:rlim')
    [~, ri, li] = decaymap(ls, xgv, ygv, zgv);
    r(:,:,:,i) = ri.data;
    l(:,:,:,i) = li.data;
    warning('on', 'pcd:mapping:rlim')
    
    % Build the local reflectivity map.
    warning('off', 'pcd:mapping:rlim')
    [~, hi, mi] = refmap(ls, xgv, ygv, zgv);
    h(:,:,:,i) = hi.data;
    m(:,:,:,i) = mi.data;
    warning('on', 'pcd:mapping:rlim')
    
    parprogress;
end
parprogress(0);

% Sum the returns and traversals up over time.
r = cumsum(r, 4);
l = cumsum(l, 4);
h = cumsum(h, 4);
m = cumsum(m, 4);

%% Save map.
save(convFile, 'dataset', 'dataFolder', 'mappingFile', 'r', 'l', ...
    'h', 'm', '-v7.3');
display(['Result written to ', lidarMapFile, '.'])
