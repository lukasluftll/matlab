% Checks whether two PCD files are equivalent, although one contains a
% point cloud in Cartesian coordinates and the other contains spherical
% coordinates.

% Load the PCD files.
pcdcart = pcdread('campus_cart.pcd');
pcdsph = pcdread('campus_sph.pcd');

% Build a point cloud from the Cartesian coordinates in the world frame.
pccart = pointCloud([pcdcart.x(:), pcdcart.y(:), pcdcart.z(:)]);

% Convert the spherical coordinates to a Cartesian point cloud in the world
% frame.
sp = trquat2tform(...
    [pcdsph.sensor_x(:), pcdsph.sensor_y(:), pcdsph.sensor_z(:)], ...
    [pcdsph.sensor_qw(:), pcdsph.sensor_qx(:), pcdsph.sensor_qy(:), ...
    pcdsph.sensor_qz(:)]);
rlim = [2,120];
ls = laserscan(sp, pcdsph.azimuth, pcdsph.elevation, pcdsph.radius, rlim);
pcsph = ls2pc(ls);

% Visualize the difference between the point clouds.
figure;
pcshowpair(pcsph, pccart);

% Plot the distances between the points of the point clouds in a histogram.
figure;
hist(pccompare(pcsph, pccart), 1000);
