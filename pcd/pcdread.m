function data = pcdread(filename)
% PCDREAD Read contents of point cloud from PCD file.
%   DATA = PCDREAD(FILENAME) reads header and data fields of the 
%   ASCII-coded PCD file FILENAME and returns them in struct DATA. 
%
%   DATA contains the following elements read from the PCD header.
%
%   'version'       PCD version
%   'width'         width of the point cloud
%   'height'        height of the point cloud
%   'viewpoint'     4x4 homogeneous transformation matrix that represents 
%                   the sensor frame defined in the PCD file
%   'count'         number of points contained in the point cloud
%
%   Moreover, DATA contains the PCD field data as matrices of size 
%   HEIGHTxWIDTHxCOUNT, with WIDTH and HEIGHT being the width and height of 
%   the point cloud and COUNT being the data element count of the 
%   respective field as specified in the PCD header.
%
%   If an accompanying DAT file is available, its contents are also read
%   and are accessible via DATA. For the naming of the respective elements 
%   of DATA, see POSREAD. 
%   In order for the DAT file to be recognized correctly, it must obey the 
%   following naming convention:
%
%      PCD file name   <filename>.pcd
%      DAT file name   <filename>_info.dat
%  
%   Example:
%      data = pcdread('office.pcd');
%      cloud = pointCloud([data.x; data.y; data.z]');
%      pcshow(cloud);
%
%   See also POSREAD, PCREAD, PTSREAD.
 
% Copyright 2016 Alexander Schaefer

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

%% Read PCD header.
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
fieldsize = cell2mat(textscan(fieldsize, '%d')); %#ok<*NASGU>

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

%% Check file contents.
% Warn if the PCD file version is older than the official entry point for
% the PCD file format.
if pcdversion < 0.7
    warning('PCD file version is below 0.7: file format not well defined.')
end

% If the number of all points does not equal the point cloud height x
% width, issue a warning.
if count ~= width * height
    warning('Inconsistent PCD header fields: POINTS ~= WIDTH * HEIGHT.')
end

% Make sure the PCD file's data type is ASCII. 
if ~strcmpi(pcdtype, 'ascii')
    error('This function only reads PCD files of type ASCII.')
end

%% Read PCD payload data.
% Read raw data.
rawdata = dlmread(filename, ' ', 11, 0);

% Convert raw data to cell array. Each cell element contains a matrix of
% size WIDTH x HEIGHT x COUNT.
celldata = cell(length(fieldcount), 1);
col = [0; cumsum(fieldcount)];
for i = 1 : length(fieldcount)
    celldata{i} = rawdata(:, col(i)+1:col(i+1));
    celldata{i} = reshape(celldata{i}, width, height, fieldcount(i));
    celldata{i} = permute(celldata{i}, [2 1 3]);
end

% If the data type of a raw data column is not double, convert it to the 
% data type given in the PCD file header.
uint = regexp(fieldtype', '[uU]');
celldata(uint) = cellfun(@uint32, celldata(uint), 'UniformOutput', false);
int = regexp(fieldtype', '[iI]');
celldata(int) = cellfun(@int32, celldata(int), 'UniformOutput', false);

% Transform the cell array into a structure array.
data = cell2struct(celldata, fieldname, 1);

%% Add header data to return structure.
data.version = pcdversion;
data.width = width;
data.height = height;
data.viewpoint = eye(4);
data.viewpoint(1:3,:) = [quat2rotm(viewpoint(4:7)'), viewpoint(1:3)];
data.count = count;

%% Read DAT file.
% Determine the name of the DAT file.
datFilename = [filename(1:end-length('.pcd')), '_info.dat'];

% Check if the DAT file exists.
if exist(datFilename, 'file') ~= 2
    return
end

% Check if the DAT file can be opened for read access.
datFid = fopen(datFilename, 'r');
if datFid == -1
    return
end
fclose(datFid);

% Read position information if available.
pos = posread(datFilename);

% Append position information to return structure.
posfield = fieldnames(pos);
for i = 1 : numel(posfield)
    data.(posfield{i}) = pos.(posfield{i});
end

end
