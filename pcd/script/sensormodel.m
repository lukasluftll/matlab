% Evaluate sensor model for BoniRob localization on leek field.

% Define parameters.
res = 0.1;
lim = [-10, 100; -10, 20];
sensorfile = 'pcd/data/sensmiddle.pcd';
rpy = [0,0,0];

% Create elevation map of field.
pcfield = pcd2pc(pcdread('pcd/data/leek.pcd'));
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));
em = elevationmap(pcfield, 0.05);

% Fill gaps in elevation map.
em = em.fillnan([5,5]);

% Read sensor measurements.
pcsens = pcd2pc(pcdread(sensorfile));
pcsens = pctransform(pcsens, ht2affine3d(eul2tform(rpy)));
psens = permute(pcsens.Location, [2,3,1]);

% Shift the measurements and compute the height difference at each point.
x = lim(1,1) : res : lim(1,2);
y = lim(2,1) : res : lim(2,2);
nx = numel(x);
ny = numel(y);
d = NaN(nx, ny);
nanfrac = NaN(nx, ny);
n = size(psens, 1);
progressbar(nx)
parfor ix = 1 : nx
    for iy = 1 : ny       
        z = mean(em-(psens+repmat([x(ix),y(iy),0], n, 1)), 'omitnan');
        dz = constrain(em-(psens+repmat([x(ix),y(iy),z], n, 1)), [0,+Inf]);
        d(ix,iy) = mean(dz, 'omitnan');
        
        nanfrac(ix,iy) = sum(isnan(dz)) / numel(dz);
    end
    progressbar
end

% Plot result.
figure('Name', 'Sensor model in middle of field')
surf(x, y, d', 'EdgeColor', 'none')
daspect([1 1 0.01])
labelaxes

% Colorize field map according to visualize the sensor model output.
c = floor(normm(d) * 63.9999) + 1;
palette = uint8(round(colormap * 255));
roi = [lim(1,1), lim(1,2)-eps, lim(2,1), lim(2,2)-eps, -Inf, +Inf];
pcfieldsel = select(pcfield, findPointsInROI(pcfield, roi));
minxy = repmat(lim(1:2,1).', pcfieldsel.Count, 1);
i = floor((pcfieldsel.Location(:,1:2) - minxy) / res) + 1;
i = sub2ind(size(c), i(:,1), i(:,2));
pcfieldsel.Color = palette(c(i),:);
figure('Name', 'Field map showing robot location probability')
pcshow(pcfieldsel, 'MarkerSize', 40)

% Show the fraction of missing correspondences.
figure('Name', 'NaN count')
surf(x, y, nanfrac', 'EdgeColor', 'none')
labelaxes

% Plot point cloud of field and the scan point cloud at the most probable
% location.
[~,imin] = min(d(:));
[xmin,ymin] = ind2sub(size(d), imin);
offset = repmat([x(xmin),y(ymin),0], n, 1);
offset(:,3) = mean(em-(psens+offset), 'omitnan');
figure('Name', 'Field')
pcshowpair(pcfield, pointCloud(psens+offset))
axis equal
labelaxes
