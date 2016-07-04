function pos = posread(filename)
% POSREAD Read position information about a laser scan.
%   POS = POSREAD(FILENAME) reads position information from the DAT file
%   specified by the string FILENAME. The return value POS is a structure 
%   with the following elements:
%
%   'time'          scalar that represents the time stamp of the laser scan
%
%   'odometry'      affine3d object that represents the odometry pose at 
%                   which the scan was captured
%
%   'gps'           structure containing the GPS coordinates of the place 
%                   where the scan was captured; its elements are:
%                   'time'            time stamp of the GPS reading
%                   'longitude'       GPS longitude
%                   'latitude'        GPS latitude
%                   'elevation'       GPS elevation
%
%   DAT file format
%   ---------------
%   The following example of a DAT file shows how the position information
%   in such a file is ordered and formatted:
%
%   Time: 1454504682.921031963
%   Odometry: -4.21839 -3.69559 0.955 0 -0 2.11513
%   GPS: 1454504682.677756915 6.96038 47.4405 0
% 
%   Example:
%      pos = posread('campus_info.dat')
%
%   See also PCDREAD, PCREAD, PTSREAD.
 
%  Copyright 2016 Alexander Schaefer

%% Validate input.
% Check the number of input arguments.
narginchk(1, 1);

% Make sure the given file name is a string.
if ~ischar(filename)
    error(message('vision:pointcloud:badFileName'));
end

% Verify that the file exists.
fid = fopen(filename, 'r');
if fid == -1
    if isempty(dir(filename))
        error(message('MATLAB:imagesci:imread:fileDoesNotExist', ...
              filename));
    else
        error(message('MATLAB:imagesci:imread:fileReadPermission', ...
              filename));
    end
end

%% Read file.
% Define the format error message.
formatError = 'Invalid input file format.';

% Read time stamp.
pos.time = 0;
line = fgetl(fid);
if ischar(line)
    time = textscan(line, '%s %f64');
    if ~strcmpi(time{1}, 'Time:')
        error(formatError);
    else
        pos.time = time{2};
    end
end

% Read odometry data.
pos.odometry = affine3d();
line = fgetl(fid);
if ischar(line)
    odometry = textscan(line, '%s %f %f %f %f %f %f');
    if ~strcmpi(odometry{1}, 'Odometry:')
        error(formatError);
    else
        odometry = [odometry{2:7}];
        translation = [odometry(1:3), 1];
        rotation = eul2rotm(odometry(4:6));
        transform(1:3, 1:3) = rotation;
        transform(4, 1:4) = translation;
        pos.odometry = affine3d(transform);
    end
end

% Read GPS data.
pos.gps = [];
line = fgetl(fid);
if ischar(line)
    gps = textscan(line, '%s %f64 %f %f %f');
    if ~strcmpi(gps{1}, 'GPS:')
        error(formatError);
    else
        pos.gps.time = gps{2};
        pos.gps.longitude = gps{3};
        pos.gps.latitude = gps{4};
        pos.gps.elevation = gps{5};
    end
end

% Close the file.
fclose(fid);

end
