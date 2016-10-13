% Evaluate elevation map sensor model for leek field dataset.

%% Define parameters.
% Step size when shifting the scan.
shiftres = 0.1;

% Minimum and maximum shifting offset in x and y direction.
shiftlim = [-10, 100; -10, 20];

% PCD file containing the map.
mapfile = 'pcd/data/leek.pcd';

% PCD file containing the scan.
sensorfile = 'pcd/data/sensmiddle.pcd';

% Orientation of the scan.
rpysens = deg2rad([0,0,0]);

%% Create elevation map of field.
% Read PCD file to point cloud.
pcfield = removeInvalidPoints(pcread(mapfile));

% Apply rotation to map.
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));

% Create the elevation map and fill the gaps in the map.
em = fillnan(elevationmap(pcfield, 0.05), [5,5]);

%% Read sensor measurements.
% Read PCD file to point cloud.
pcsens = removeInvalidPoints(pcdread(sensorfile));

% Rotate the scan.
pcsens = pctransform(pcsens, ht2affine3d(eul2tform(rpysens)));

%% Shift scan.
% Define the offset vectors.
x = shiftlim(1,1) : shiftres : shiftlim(1,2);
y = shiftlim(2,1) : shiftres : shiftlim(2,2);
nx = numel(x);
ny = numel(y);

% Initialize the matrix that stores the mean z difference.
d = NaN(nx, ny);

% Initialize the matrix that stores the fraction of NaN differences.
nanfrac = NaN(nx, ny);

% Vary the position of the scan and compute the z difference for
% each offset.
progressbar(nx)
parfor ix = 1 : nx
    for iy = 1 : ny
        % Extract the Mx3 vector of point coordinates from the scan.
        psens = reshape(pcsens.Location(:),pcsens.Count,3,1); %#ok<*PFBNS>
        
        % Shift the scan horizontally.
        psens = psens + repmat([x(ix),y(iy),0], size(psens,1), 1);
        
        % Adjust the z coordinate of the scan origin.
        ifin = all(isfinite(psens(:,1:2)), 2);
        psens = psens - mean(em.diff(psens), 'omitnan');
        
        % Compute the difference in z for all points of the scan.
        % Allow only positive values, as negative differences mean the scan
        % hit an object and the corresponding point is correct.
        dz = constrain(em.diff(psens), [0,+Inf]);
        
        % Compute the fraction of unmatched points.
        nanfrac(ix,iy) = sum(isnan(dz)) / numel(dz);
        
        % Compute the mean difference along z between the elevation map and 
        % the scan.
        if nanfrac(ix,iy) < 0.1
            d(ix,iy) = mean(dz, 'omitnan');
        end        
    end
    
    % Advance the progress bar.
    progressbar
end

%% Plot results.
% Plot the mean z distance for all offsets.
fig = figure('Name', 'Mean distance in z');
surf(x, y, d', 'EdgeColor', 'none')
daspect([1,1,0.01])
labelaxes
