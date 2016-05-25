function position = posread(filename)
% posread Read a file containing position information about a laser scan.
%   position = posread(filename) reads position information from the 
%   DAT file specified by the string filename. The return value 
%   position is a structure with the following elements:
%
%   time            scalar that represents the time stamp of the laser scan
%
%   odometry        affine3d object that represents the odometry pose at 
%                   which the scan was captured.
%
%   gps             structure containing the GPS coordinates of the place 
%                   where the scan was captured; the structure's elements
%                   are:
%                   time            time stamp of the GPS reading
%                   longitude       GPS	longitude
%                   latitude        GPS latitude
%                   elevation       GPS elevation
%                   
%
%   DAT file format
%   ---------------
%   The following example of a DAT file shows how the position information
%   in such a file is formatted:
%
%   Time: 1454504682.921031963
%   Odometry: -4.21839 -3.69559 0.955 0 -0 2.11513
%   GPS: 1454504682.677756915 6.96038 47.4405 0
% 
%   Example : read a position from a DAT file
%   -----------------------------------------
%   position = posread('terrain.dat');
%
%   See also pcread, ptsread, pcdread, pointCloud.
 
%  Copyright 2016 Alexander Schaefer

%% Validate the input.
% Make sure the given file name is a string.
if ~ischar(filename)
    error(message('vision:pointcloud:badFileName'));
end

% Verify that the file exists.
fid = fopen(filename, 'r');
if (fid == -1)
    if isempty(dir(filename))
        error(message('MATLAB:imagesci:imread:fileDoesNotExist', ...
              filename));
    else
        error(message('MATLAB:imagesci:imread:fileReadPermission', ...
              filename));
    end
end

%% Read the file entries.
% Define the format error message.
formatError = 'Invalid input file format.';

% Read the time stamp.
position.time = 0;
time = textscan(fgetl(fid), '%s %f64');
if (~strcmpi(time{1}, 'Time:'))
    error(formatError);
else
    position.time = time{2};
end

% Read the odometry data.
position.odometry = affine3d();
odometry = textscan(fgetl(fid), '%s %f %f %f %f %f %f');
if (~strcmpi(odometry{1}, 'Odometry:'))
    error(formatError);
else
    odometry = [odometry{2:7}];
    translation = [odometry(1:3), 1];
    rotation = eul2rotm(odometry(4:6));
    transform = zeros(4);
    transform(1:3, 1:3) = rotation;
    transform(4, :) = translation;
    position.odometry = affine3d(transform);
end

% Read the GPS data.
position.gps = [];
gps = textscan(fgetl(fid), '%s %f64 %f %f %f');
if (~strcmpi(gps{1}, 'GPS:'))
    error(formatError);
else
    position.gps.time = gps{2};
    position.gps.longitude = gps{3};
    position.gps.latitude = gps{4};
    position.gps.elevation = gps{5};
end

end
