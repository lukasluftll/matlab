function ls = lsread(file, rlim)
% LSREAD Read PCD file to laserscan object.
%   LS = LSREAD(FILE) uses the information contained in the PCD file FILE 
%   to create a laserscan object LS.
%
%   FILE is the full path to a PCD file.
%
%   LS =  LSREAD(FILE, RLIM) additionally passes the ordered 2-element 
%   vector RLIM to the constructor of laserscan to specify the lidar sensor
%   reading range.
%
%   Example:
%      ls = lsread('pcd/data/campus.pcd', [2,120])
%
%   See also PCDREAD, PCDREAD.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check number of input arguments.
narginchk(1, 2)

% Check the input arguments.
validateattributes(file, {'char'}, {'row', 'nonempty'}, '', 'FILE')
if nargin > 1
    validateattributes(rlim, {'numeric'}, {'numel', 2}, '', 'RLIM')
end

%% Read PCD file.
pcd = pcdread(file);

% Map all infinite radius values to NaN.
pcd.radius(~isfinite(pcd.radius)) = NaN;

%% Construct laserscan object.
% Build a 3D matrix that contains the homogeneous coordinate
% transformations indicating the pose of the sensor for each ray of the
% scan.
sp = trquat2tform([pcd.sensor_x(:), pcd.sensor_y(:), pcd.sensor_z(:)], ...
    [pcd.sensor_qw(:),pcd.sensor_qx(:),pcd.sensor_qy(:),pcd.sensor_qz(:)]);

% Create the laserscan object.
if nargin < 2
    ls = laserscan(sp, pcd.azimuth(:), pcd.elevation(:), pcd.radius(:));
else
    ls = laserscan(sp,pcd.azimuth(:),pcd.elevation(:),pcd.radius(:),rlim);
end

end
