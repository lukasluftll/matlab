function data = pcdread(filename)
% PCDREAD Read contents of point cloud from PCD file.
%   DATA = PCDREAD(FILENAME) reads all data fields of the PCD file defined 
%   by FILENAME and returns them in a struct DATA.
%   If an accompanying DAT file is available, its contents are also added
%   to DATA.
%
%   DAT file naming convention
%   --------------------------
%   The DAT file must obey the following naming convention:
% 
%   PCD file name   <filename>.pcd
%   DAT file name   <filename>_info.dat
%
%   Recognized field names
%   ----------------------
%   PCDREAD recognizes the following data field names in the PCD file 
%   header and returns the corresponding data types:
%
%   'viewpoint'     affine3d object that represents the sensor frame 
%                   defined in the PCD file
%
%   'odometry'      affine3d object that represents the odometry pose 
%                   defined in the DAT file
%
%   'gps'           structure array that represents the GPS coordinates 
%                   of the DAT file; it contains the following elements:
%                   'time'            time stamp of the GPS reading
%                   'longitude'       GPS	longitude
%                   'latitude'        GPS latitude
%                   'elevation'       GPS elevation
% 
%   Organized vs. unorganized point clouds
%   --------------------------------------
%   If the PCD file header specifies a height value of 1, the point cloud
%   is unorganized. In this case, the values of the data fields are stored 
%   in matrices of size 1xWIDTHxCOUNT, with WIDTH being the width of the 
%   point cloud and COUNT being the data element count from the PCD file 
%   header. 
%
%   If the specified height value is greater than 1, the point cloud is
%   organized and the values are stored in matrices of size
%   HEIGHTxWIDTHxCOUNT.
% 
%   PCD file format
%   ---------------
%   PCD files occur in different formats: ASCII and binary. This function
%   only reads ASCII files.
%  
%   Example:
%      data = pcdread('office.pcd');
%      cloud = pointCloud([data.x; data.y; data.z]');
%      pcshow(cloud);
%
%   See also PCREAD, PTSREAD, POINTCLOUD, PCWRITE, PCSHOW.
 
% Copyright 2016 Alexander Schaefer

%% Validate PCD file name.
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

%% Read PCD header entries.
% Skip the opening line.
fgetl(fid);

% PCD file version.
pcdversion = fgetl(fid);
pcdversion = pcdversion(length('VERSION')+1:end);
pcdversion = cell2mat(textscan(pcdversion, '%f'));

% Fields contained in the PCD file.
fieldname = fgetl(fid);
fieldname = fieldname(length('FIELDS')+1:end);
fieldname = textscan(fieldname, '%s');
fieldname = fieldname{:};

% Size of each field in bytes.
fieldsize = fgetl(fid);
fieldsize = fieldsize(length('SIZE')+1:end);
fieldsize = cell2mat(textscan(fieldsize, '%d'));

% Data type of each field.
fieldtype = fgetl(fid);
fieldtype = fieldtype(length('TYPE')+1:end);
fieldtype = cell2mat(textscan(fieldtype, '%c'));

% Number of elements per field per point.
fieldcount = fgetl(fid);
fieldcount = fieldcount(length('COUNT')+1:end);
fieldcount = cell2mat(textscan(fieldcount, '%d'));

% Point cloud width.
width = fgetl(fid);
width = width(length('WIDTH')+1:end);
width = cell2mat(textscan(width, '%d'));

% Point cloud height.
height = fgetl(fid);
height = height(length('HEIGHT')+1:end);
height = cell2mat(textscan(height, '%d'));

% Acquisition viewpoint of the points in the dataset.
viewpoint = fgetl(fid);
viewpoint = viewpoint(length('VIEWPOINT')+1:end);
viewpoint = cell2mat(textscan(viewpoint, '%f'));

% Total number of points.
count = fgetl(fid);
count = count(length('POINTS')+1:end);
count = cell2mat(textscan(count, '%d'));

% Data type of the PCD file.
pcdtype = fgetl(fid);
pcdtype = pcdtype(length('DATA')+1:end);
pcdtype = textscan(pcdtype, '%s');
pcdtype = pcdtype{:};

fclose(fid);

%% Check file and read payload data.
% Warn if the PCD file version is older than the official entry point for
% the PCD file format.
if pcdversion < 0.7
    warning('PCD file version is below 0.7: file format not well defined.')
end

% Make sure the PCD file's data type is ASCII. 
if ~strcmpi(pcdtype, 'ascii')
    error('This function only reads PCD files of type ''ascii''.')
end

% If the number of all points does not equal the point cloud height x
% width, issue a warning.
if count ~= width * height
    warning('Inconsistent PCD header fields: POINTS ~= WIDTH * HEIGHT.')
end

% Read raw data.
rawdata = dlmread(filename, ' ', 11, 0);

%% Read payload data.
% Convert raw data to cell array. Each cell element contains a matrix of
% size WIDTH x HEIGHT x FIELDCOUNT.
data = cell(1, length(fieldcount));
col = [0; cumsum(fieldcount)];
for i = 1 : length(fieldcount)
    data{i} = rawdata(:, col(i)+1:col(i+1));
    data{i} = reshape(data{i}, width, height, fieldcount(i));
    data{i} = permute(data{i}, [2 1 3]);
end

% Convert data to the data types specified in the PCD file header.
uint = regexp(fieldtype', '[uU]');
data(uint) = cellfun(@uint32, data(uint), 'UniformOutput', false);
int = regexp(fieldtype', '[iI]');
data(int) = cellfun(@int32, data(int), 'UniformOutput', false);

% Transform the cell array into a structure array.
data = cell2struct(data', fieldname, 1);

%% Append viewpoint to return structure.
viewpointMatrix = quat2rotm(viewpoint(4:7)');
viewpointMatrix(:, 4) = viewpoint(1:3)';
viewpointMatrix(4, :) = [0, 0, 0, 1];
data.viewpoint = affine3d(viewpointMatrix');

%% Read DAT file and append information to return structure.
% Create the DAT file name following naming convention.
datFilename = [filename(1:end-length('.pcd')), '_info.dat'];

% Read position information.
fid = fopen(datFilename, 'r');
if fid ~= -1
    pos = posread(datFilename);
    
    % Append position information to return structure.
    posfield = fieldnames(pos);
    for i = 1 : numel(posfield)
        data.(posfield{i}) = pos.(posfield{i});
    end
end

end
