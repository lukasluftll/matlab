% Voxelizes the volume occupied by a Lidar scan and computes for each
% voxel the reflection probability.

% Define the voxel map resolution.
res = 1.0;

% Read the point cloud.
cloud = pcdread('data/office.pcd');
pointsCart = cloud.pointCloud.Location;

% Compute the spherical coordinates of the points w.r.t. the sensor frame.
beamsCart = pointsCart ...
    - repmat(cloud.viewpoint.T(4,1:3), size(pointsCart, 1), 1);
[azimuth, elevation, r] = cart2sph(...
    beamsCart(:,1), beamsCart(:,2), beamsCart(:,3));
beamsSph = [azimuth, elevation, r];

% Compute the volume occupied by the Lidar scan.
limits = [cloud.pointCloud.XLimits; ...
    cloud.pointCloud.YLimits; ...
    cloud.pointCloud.ZLimits];
limits = [floor(limits(:,1)), ceil(limits(:,2))];    

% Create the figure where to plot the occupancy map.
fig = figure('Name', 'Occupancy map', 'NumberTitle', 'Off');

% Iterate over all voxels.
for (x = limits(1,1) : res : limits(1,2))
    for (y = limits(2,1) : res : limits(2,2))
        for (z = limits(3,1) : res : limits(3,2))
            % Compute the limits of the voxel.
            roi = [x, x+res; y, y+res; z, z+res];
            
            % Find all points living in the voxel.
            roiIndices = findPointsInROI(cloud.pointCloud, roi);
            roiCartPoints = pointsCart(roiIndices,:);
            hits = length(roiIndices);
            
            % Compute the coordinates of all 8 voxel corners.
            voxelCorners = [x, y, z];
            voxelCorners = [voxelCorners; ...
                voxelCorners + repmat([res, 0, 0], 1, 1)];
            voxelCorners = [voxelCorners; ...
                voxelCorners + repmat([0, res, 0], 2, 1)];
            voxelCorners = [voxelCorners; ...
                voxelCorners + repmat([0, 0, res], 4, 1)];
            
            % Use the voxel corners to compute the minimum and maximum 
            % azimuth and elevation angles of beams that traverse the 
            % voxel.
            [azimuth, elevation, r] = cart2sph(...
                voxelCorners(:,1), voxelCorners(:,2), voxelCorners(:,3));
            
            % Make sure the azimuth difference for the voxel corners is
            % correctly computed if the voxel touches the negative x-axis;
            % there, the azimuth angle changes from +pi to -pi.
            if (min(voxelCorners(:,1)) < 0 && max(voxelCorners(:,2)) == 0)
                azimuth(azimuth == pi) = -pi;
            end
            
            % Delete all beams that cannot dip in the voxel volume.
            travbeams = beamsSph;
            travbeams(travbeams(:,1) < min(azimuth), :) = [];
            travbeams(travbeams(:,1) > max(azimuth), :) = [];
            travbeams(travbeams(:,2) < min(elevation), :) = [];
            travbeams(travbeams(:,2) > max(elevation), :) = [];
            travbeams(travbeams(:,3) < min(r), :) = [];
            
            % Find out how many beams dipped into the voxel volume.
            dips = 0;
            for (b = 1 : size(travbeams, 1))
                % Compute the intersection of the beam and the voxel.
                origin = zeros(1, 3);
                direction = sph2cart(travbeams(b,:));
                vmin = roi(:,1);
                vmax = roi(:,2);
                [flag, tmin] = rayBoxIntersection(...
                    origin, direction, vmin, vmax);
                
                % Check whether the beam dipped into the voxel volume.
                if (tmin < 1)
                    dips = dips + 1;
                end
            end
            
            % Plot the cube.
            if (dips < 1)
                cube([x + res/2, y + res/2, z + res/2], res);
            end            
        end
    end
end
