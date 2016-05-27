% Voxelizes the volume occupied by a Lidar scan and computes for each
% voxel the reflection probability.
% The reflection probability is computed under the assumption that the 
% points in each voxel are Gauss-distributed.

% Define the voxel map resolution.
res = 1;

% Read the point cloud.
cloud = pcdread('data/campus.pcd');
points = cloud.pointCloud.Location;
pointcount = numel(points(:,:,1));

% Compute the spherical coordinates of the points w.r.t. the sensor frame.
beams = pc2sph(points);

% Compute the volume occupied by the Lidar scan.
limits = [cloud.pointCloud.XLimits; ...
    cloud.pointCloud.YLimits; ...
    cloud.pointCloud.ZLimits];

% Make sure the limits are multiples of the resolution so the voxels
% are alinged along the planes spanned by the axis pairs.
limits = [floor(limits(:,1) / res) * res, ceil(limits(:,2) / res) * res];

% Compute the diameter of the point cloud.
diameter = norm(limits(:,2) - limits(:,1));

% Create the figure where to plot the occupancy map.
fig = figure('Name', 'Occupancy map', 'NumberTitle', 'Off');
grid on;
axis equal;

% Plot the point cloud.
pcshow(cloud.pointCloud, 'MarkerSize', 80);
drawnow;

% Iterate over all voxels and for each one compute the reflection 
% probability.
voxelcount = (limits(:,2) - limits(:,1))' / res;
mu = zeros([3, voxelcount]);
cov = zeros([3, 3, voxelcount]);
alpha = 0.5 * ones(voxelcount);
ix = 1;
for x = limits(1,1) : res : limits(1,2)
    iy = 1;
    for y = limits(2,1) : res : limits(2,2)
        iz = 1;
        for z = limits(3,1) : res : limits(3,2)
            % Compute the limits of the voxel.
            voxel = [x, x+res; y, y+res; z, z+res];
            
            % Find all points that reside inside the voxel.
            roiIndices = findPointsInROI(cloud.pointCloud, voxel);
            hits = numel(roiIndices);
            roiPoints = [points(roiIndices) + 0*pointcount, ...
                points(roiIndices + 1*pointcount), ...
                points(roiIndices + 2*pointcount)];
            
            % Compute the coordinates of all 8 voxel vertices.
            vertices = [x, y, z];    %#ok<*AGROW>
            vertices = [vertices; vertices + repmat([res, 0, 0], 1, 1)]; 
            vertices = [vertices; vertices + repmat([0, res, 0], 2, 1)];
            vertices = [vertices; vertices + repmat([0, 0, res], 4, 1)];
            
            % Use the vertices to compute the minimum and maximum 
            % azimuth and elevation angles of the beams that penetrate
            % the voxel.
            [azimuth, elevation, radius] = cart2sph(...
                vertices(:,1), vertices(:,2), vertices(:,3));
            
            % Make sure the azimuth difference for the voxel corners is
            % correctly computed if the voxel touches the negative x-axis;
            % at the negative x-axis, the azimuth angle changes 
            % from +pi to -pi.
            if (min(vertices(:,1)) < 0 && max(vertices(:,2)) == 0)
                azimuth(azimuth == pi) = -pi;
            end
            
            % Find the logical subscripts of all beams that can
            % permeate the voxel volume.
            keep = beams(:,:,1) >= min(azimuth) ...
                & beams(:,:,1) <= max(azimuth) ...
                & beams(:,:,2) >= min(elevation) ...
                & beams(:,:,2) <= max(elevation) ...
                & beams(:,:,3) >= min(radius);
            
            % Consider only those beams in the following.
            beamsAzimuth = beams(:,:,1);
            beamsElevation = beams(:,:,2);
            beamsRadius = beams(:,:,3);
            permbeams = cat(3, beamsAzimuth(keep), ...
                beamsElevation(keep), beamsRadius(keep));
            
            % Find out how many beams really permeate the voxel volume.
            perms = 0;
            for row = 1 : size(permbeams, 1)
                for col = 1 : size(permbeams, 2)
                    % NaN beams are infinitely long; set their radius to 
                    % the diameter of the point cloud.
                    if isnan(permbeams(row,col,3))
                        permbeams(row,col,3) = diameter;
                    end

                    % Compute the intersection of the beam and the voxel.
                    origin = zeros(1, 3);
                    [dirx, diry, dirz] = sph2cart(permbeams(row,col,1), ...
                        permbeams(row,col,2), permbeams(row,col,3));
                    vmin = voxel(:,1);
                    vmax = voxel(:,2);
                    [flag, tmin] = rayBoxIntersection(...
                        origin, [dirx, diry, dirz], vmin, vmax);

                    % Check whether the beam permeates the voxel volume.
                    if (tmin < 1)
                        perms = perms + 1;
                    end
                end
            end
            
            % Plot the voxel.
            if (perms > 0)
                center = [x + res/2, y + res/2, z + res/2];
                
                % If at least one beam permeated the voxel, plot it red.
                % The transparency indicates the probability of a
                % reflection.
                % If the voxel did never reflect a beam, color it gray.
                 if (hits > 0)
                     alpha(ix,iy,iz) = hits / perms;
                     cube(center, res, 'FaceColor', 'red', ...
                         'FaceAlpha', alpha(ix,iy,iz)*0.8 + 0.2, ...
                         'LineStyle', 'none');
                 else
                    cube(center, res, 'FaceColor', [0, 0, 0], ...
                        'FaceAlpha', 0.05, 'LineStyle', 'none');
                 end
            end
            
            % Calculate the properties of the Gaussian distribution of
            % the points.
            mu(:,ix,iy,iz) = mean(roiPoints, 1);
            cov(:,:,ix,iy,iz) = roiPoints' * roiPoints;
            
            iz = iz + 1;
        end
        iy = iy + 1;
    end
    ix = ix + 1;
end

% Create a figure in which to plot the reflection probabilities.
figure('Name', 'Reflection probability', 'NumberTitle', 'Off');

% Compute the coordinates of the centers of all voxels.
limits = limits + res/2;
[centerX, centerY, centerZ] = meshgrid(limits(1,1) : res : limits(1,2), ...
    limits(2,1) : res : limits(2,2), ...
    limits(3,1) : res : limits(3,2));
center = [reshape(centerX, numel(centerX), 1), ...
    reshape(centerY, numel(centerY), 1), ...
    reshape(centerZ, numel(centerZ), 1)];

% Compute the reflection probability for each voxel.
% for ix = 1 : numel(limits(1,1) : res 
%     for iy = 1 : size(centerY, 
% mvnpdf(center, mu, cov)
