% Draws two figures: one of a single Velodyne VLP-16 scan and another 
% of a bundle of scans.

% Define parameters.
res = 0.3;
palette = 'gray';
bg = [0.2, 0.2, 0.2];
markerSize = 70;

% Read a subset of the scans recorded with BoniRob on campus.
datadir = '~/ros/datasets/bonirob/campus_pcd/';
files = dir([datadir, 'velodyne_freiburg-04*.pcd']);

% Read the first point cloud.
campus = pcdread([datadir, files(1).name]);
campus = pointCloud(cat(3, campus.x, cat(3, campus.y, campus.z)));

% Display the single scan.
fig = figure('Name', 'Single scan', 'NumberTitle', 'Off');
singleax = pcshow(campus, 'MarkerSize', markerSize);
whitebg(fig, bg);

% Merge the scans to form a bundle.
for (i = 2 : min(numel(files), 21))
    pc = pcdread([datadir, files(i).name]);
    pc = pointCloud(cat(3, pc.x, cat(3, pc.y, pc.z)));
    campus = pcmerge(campus, pc, res);
end

% Display the scan bundle.
fig = figure('Name', 'Scan bundle', 'NumberTitle', 'Off');
bundleax = pcshow(campus, 'MarkerSize', markerSize);
whitebg(fig, bg);

% Link the cameras of the two figures together.
cameralink = linkprop([singleax, bundleax], {'CameraUpVector', ...
    'CameraPosition', 'CameraTarget', 'CameraViewAngle'});
setappdata(gcf, 'StoreTheLink', cameralink);
