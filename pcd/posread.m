function pos = posread(filename)
% POSREAD Read position information from file.
%   POS = POSREAD(FILENAME) reads position information from the DAT file
%   specified by the string FILENAME. The return value POS is a structure 
%   that contains the following elements.
%
%   'time'          scalar that represents the time stamp of the laser scan
%
%   'odometry'      4x4 homogeneous transformation matrix that represents 
%                   the odometry pose at which the scan was captured
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
%   Odometry: -9.23 -62.10 -0.11 0.787 -0.006 -0.025 0.617
%   GPS: 1454504682.677756915 6.96038 47.4405 0
%
%   Odometry information is given as a list of Cartesian position and 
%   rotation as quaternion: x, y, z, qw, qx, qy, qz.
%   GPS information is given as a list of time, longitude, latitude, 
%   elevation.
%
%   If a line is missing in the DAT file, the corresponding element of POS
%   is empty.
% 
%   Example:
%      pos = posread('campus_info.dat')
%
%   See also PCDREAD, PCREAD, PTSREAD.
 
%  Copyright 2016 Alexander Schaefer

%% Validate input.
% Check the number of input arguments.
narginchk(1, 1)

% Make sure the given file name is a string.
if ~ischar(filename)
    error('FILENAME must be a string.')
end

% Verify the file exists.
if exist(filename, 'file') ~= 2
    error(['File ', filename, ' does not exist.'])
end

%% Open file.
% Try to open the file.
fid = fopen(filename, 'r');
if fid == -1
    error(['Cannot read ', filename, '.'])
end

% Make sure the file is closed upon function termination.
cleaner = onCleanup(@() fclose(fid));

%% Read file content.
% Define the format error message.
formatError = 'Invalid input file format.';

% Read time stamp.
pos.time = [];
line = fgetl(fid);
if ischar(line)
    time = textscan(line, '%s %f64');
    if strcmpi(time{1}, 'Time:')
        pos.time = time{2};
    else
        error(formatError);        
    end
end

% Read odometry data.
pos.odometry = [];
line = fgetl(fid);
if ischar(line)
    odometry = textscan(line, '%s %f %f %f %f %f %f %f');
    if strcmpi(odometry{1}, 'Odometry:')
        pos.odometry = quat2tform([odometry{5:8}]);
        pos.odometry(1:3,4) = [odometry{2:4}]';
    else
        error(formatError);
    end
end

% Read GPS data.
pos.gps = [];
line = fgetl(fid);
if ischar(line)
    gps = textscan(line, '%s %f64 %f %f %f');
    if strcmpi(gps{1}, 'GPS:')
        pos.gps.time = gps{2};
        pos.gps.longitude = gps{3};
        pos.gps.latitude = gps{4};
        pos.gps.elevation = gps{5};
    else
        error(formatError);
    end
end

end
