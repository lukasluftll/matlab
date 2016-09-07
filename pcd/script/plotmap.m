% Plot map with robot trajectory.

%% Fetch parameters.
lidarparams

%% Load point cloud.
load(pcMapFile, 'pcMap')

%% Read robot trajectory.
rt = [];
parprogress(numel(mappingFile));
for i = 1 : numel(mappingFile)
    ls = lsread([dataFolder, '/', mappingFile(i).name], rlim);
    rt = unique([rt; unique(tform2trvec(ls.sp), 'rows')], 'rows'); 
    parprogress;
end
parprogress(0);

%% Plot map.
colormap('hot')
pcshow(pcMap)
hold on
plot3(rt(:,1), rt(:,2), rt(:,3)+1, '.')
labelaxes('m')
hold off
