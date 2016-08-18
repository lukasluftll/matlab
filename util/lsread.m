function ls = lsread(file, rlim)

narginchk(1, 2)

if ~ischar(file)
    error('FILE must be a string.')
end

% Read the PCD file.
pcd = pcdread(file);

sp = trquat2tform([pcd.sensor_x(:), pcd.sensor_y(:), pcd.sensor_z(:)], ...
    [pcd.sensor_qw(:),pcd.sensor_qx(:),pcd.sensor_qy(:),pcd.sensor_qz(:)]);

if nargin < 2
    ls = laserscan(sp, pcd.azimuth(:), pcd.elevation(:), pcd.radius(:));
else
    ls = laserscan(sp,pcd.azimuth(:),pcd.elevation(:),pcd.radius(:),rlim);
end

end
