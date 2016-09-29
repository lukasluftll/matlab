% Evaluate sensor model for BoniRob localization on leek field.

% Define parameters.
res = 0.1;
lim = [-10, 100; -10, 20];

% Create elevation map of field.
pcfield = pcd2pc(pcdread('pcd/data/leek.pcd'));
pcfield = pctransform(pcfield, ht2affine3d(eul2tform([pi,0,0])));
em = elevationmap(pcfield, 0.05);

% Fill gaps in elevation map.
em = em.fillnan([5,5]);

% Read sensor measurements.
pcsens = pcd2pc(pcdread('pcd/data/sensmiddle.pcd'));
psens = permute(pcsens.Location, [2,3,1]);

% Shift the measurements and evaluate the height difference at each point.
x = lim(1,1) : res : lim(1,2);
y = lim(2,1) : res : lim(2,2);
nx = numel(x);
ny = numel(y);
d = NaN(nx, ny);
n = pcsens.Count;
progressbar(nx)
parfor ix = 1 : nx
    for iy = 1 : ny       
        offset = repmat([x(ix),y(iy),0], n, 1); %#ok<PFBNS>
        offset(:,3) = mean(em-(psens+offset), 'omitnan');
        
        dz = constrain(em-(psens+offset), [0,+Inf]);
        d(ix,iy) = mean(abs(dz), 'omitnan');
    end
    progressbar
end

% Plot result.
figure('Name', 'Sensor model in middle of field')
surf(y, x, d, 'EdgeColor', 'none')
labelaxes

% Colorize field map according to visualize the sensor model output.
c = floor(normm(d) * 63.9999) + 1;
palette = uint8(round(colormap('hot') * 255));
roi = [lim(1,1), lim(1,2)-eps, lim(2,1), lim(2,2)-eps, -Inf, +Inf];
pcfield = select(pcfield, findPointsInROI(pcfield, roi));
minxy = repmat(lim(1:2,1).',pcfield.Count,1);
i = floor((pcfield.Location(:,1:2) - minxy) / res) + 1;
i = sub2ind(size(c), i(:,1), i(:,2));
pcfield.Color = palette(c(i),:);
figure('Name', 'Field map showing robot location probability')
pcshow(pcfield, 'MarkerSize', 40)
labelaxes
