function outcloud = pcdread(filename)
% pcdread Read a 3D point cloud from PCD file.
%   outcloud = pcdread(filename) reads a point cloud from the 
%   PCD file specified by the string filename.
%   If a accompanying DAT file is available, this function also reads 
%   the DAT file's contents.
%
%   In order for the DAT file to be read, both files must stay with the 
%   following naming convention:
% 
%   PCD file name   <name>.pcd
%   DAT file name   <name>_info.dat
%
%   The return value outcloud is a structure with the following elements:
%
%   pointCloud      pointCloud object read from the PCD file
%
%   intensity       matrix that contains the intensity values stored in the 
%                   field 'intensity' in the PCD file
%
%   viewpoint       affine3d object that describes the sensor coordinate
%                   frame stores in the PCD file
%  
%   time            time stamp stored in the DAT file
%
%   odometry        affine3d object that describes the odometry pose stored
%                   in the DAT file
%
%   gps             structure that stores the GPS coordinate read from the
%                   DAT file; its elements are:
%                   time            time stamp of the GPS reading
%                   longitude       GPS	longitude
%                   latitude        GPS latitude
%                   elevation       GPS elevation
%
%   PCD file format
%   ---------------
%   PCD files occur in different formats: ASCII and binary. This function
%   only reads ASCII files.
%
%   Relevant PCD file contents
%   --------------------------
%   PCD files as defined in the Point Cloud Library can contain numerous 
%   data entries. pcdread loads the following fields, if given: 
%   x, y, z                       positions of the points
%   r, g, b                       colors of the points
%   normal_x, normal_y, normal_z  normals of the points
%   intensity                     remission values of the points
%
%   The data entries can be accessed as follows:
%   x, y, z                       outcloud.pointCloud.Location
%   r, g, b                       outcloud.pointCloud.Color
%   normal_x, normal_y, normal_z  outcloud.pointCloud.Normal
%   intensity                     outcloud.intensity
%
%   Organized vs. unorganized point clouds
%   --------------------------------------
%   If the PCD file header specifies a height value of 1, the point cloud
%   is unorganized. In this case, the points' positions, colors, and 
%   normals are stored in matrices of size [height * width, 3]. 
%   intensity is a matrix of size [height * width, 1].
%
%   If the specified height value is greater than 1, the point cloud is
%   organized and the values are stored in matrices of size 
%   [height, width, 3] and [height, width, 1], respectively.
% 
%   Example : read a point cloud from a PCD file
%   --------------------------------------------
%   outcloud = pcdread('terrain.pcd');
%   pcshow(outcloud.pointCloud);
%
%   See also pcread, ptsread, pointCloud, pcwrite, pcshow.
 
%  Copyright 2016 Alexander Schaefer

%% Validate PCD file name.
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

%% Read PCD header entries.
% Skip the opening line.
fgetl(fid);

% PCD file version.
pcdVersionString = fgetl(fid);
pcdVersionString = pcdVersionString(length('VERSION')+1 : end);
pcdVersion = textscan(pcdVersionString, '%f');
pcdVersion = pcdVersion{1};

% Fields contained in the PCD file.
fieldsString = fgetl(fid);
fieldsString = fieldsString(length('FIELDS')+1 : end);
fields = textscan(fieldsString, '%s');
fields = fields{1};

% Size of each field in bytes.
fieldSizeString = fgetl(fid);
fieldSizeString = fieldSizeString(length('SIZE')+1 : end);
fieldSize = textscan(fieldSizeString, '%d');
fieldSize = fieldSize{1};

% Data type of each field.
fieldTypeString = fgetl(fid);
fieldTypeString = fieldTypeString(length('TYPE')+1 : end);
fieldType = textscan(fieldTypeString, '%c');
fieldType = fieldType{1};

% Number of elements per dimension.
countString = fgetl(fid);
countString = countString(length('COUNT')+1 : end);
count = textscan(countString, '%d');
count = count{1};

% Point cloud width.
widthString = fgetl(fid);
widthString = widthString(length('WIDTH')+1 : end);
width = textscan(widthString, '%d');
width = width{1};

% Point cloud height.
heightString = fgetl(fid);
heightString = heightString(length('HEIGHT')+1 : end);
height = textscan(heightString, '%d');
height = height{1};

% Acquisition viewpoint for the points in the dataset.
viewpointString = fgetl(fid);
viewpointString = viewpointString(length('VIEWPOINT')+1 : end);
viewpoint = textscan(viewpointString, '%f');
viewpoint = viewpoint{1};

% Total number of points.
pointsString = fgetl(fid);
pointsString = pointsString(length('POINTS')+1 : end);
points = textscan(pointsString, '%d');
points = points{1};

% Data type of the PCD file.
dataString = fgetl(fid);
dataString = dataString(length('DATA')+1 : end);
data = textscan(dataString, '%s');
data = data{1};

fclose(fid);

%% Check file and read payload data.
% Warn if the PCD file version is older than the official entry point for
% the PCD file format.
if (pcdVersion < 0.7)
    warning('PCD file version is below 0.7 and thus not well defined.');
end

% Make sure the PCD file's data type is ASCII. 
if (~strcmpi(data, 'ascii'))
    error('This function can only read PCD files of type ASCII.');
end

% Read the data.
data = dlmread(filename, ' ', 11, 0);

%% Read sensor viewpoint.
% Convert the quaternion to a rotation matrix.
viewpointMatrix = quat2rotm(viewpoint(4:7)');

% Add the translation vector.
viewpointMatrix(:, 4) = viewpoint(1:3)';
viewpointMatrix(4, :) = [0, 0, 0, 1];

% Construct the affine3d object.
viewpoint = affine3d(viewpointMatrix');

%% Read locations.
% Check if the data contains the point coordinates.
if (isempty(strmatch('x', fields)) ...
        || isempty(strmatch('y', fields)) ...
        || isempty(strmatch('z', fields)))
    error('PCD file does not contain fields named x, y, z.');
end

% Get the coordinates of the points.
x = data(:, strmatch('x', fields)); %#ok<*MATCH2>
y = data(:, strmatch('y', fields));
z = data(:, strmatch('z', fields));

% In case the point cloud is ordered, save the point coordinates to a 
% 3D matrix.
if (height > 1)
    x = reshape(x, height, width); 
    y = reshape(y, height, width);
    z = reshape(z, height, width);
    location = cat(3, x, y, z);
else
    location = [x, y, z];
end

%% Read intensities.
% Check if the file contains intensity values.
intensity = [];
if (~isempty(strmatch('intensity', fields)))
    % Read the intensity values.
    intensity = data(:, strmatch('intensity', fields));
    
    % Organize the values, if the point cloud is organized.
    if (height > 1)
        intensity = reshape(intensity, width, height)';
    end
end

%% Read colors.
% Check if the file contains color values.
color = [];
if (~isempty(strmatch('r', fields)) ...
        && ~isempty(strmatch('g', fields)) ...
        && ~isempty(strmatch('b', fields)))
    % Read the color values.
    r = uint8(data(:, strmatch('r', fields)));
    g = uint8(data(:, strmatch('g', fields)));
    b = uint8(data(:, strmatch('b', fields)));
    
    % Organize the values, if the point cloud is organized.
    if (height > 1)
        r = reshape(r, height, width);
        g = reshape(g, height, width);
        b = reshape(b, height, width);
        color = cat(3, r, g, b);
    else
        color = [r, g, b];
    end
end

%% Read normals.
% Check if the file contains normals.
normal = [];
if (~isempty(strmatch('normal_x', fields)) ...
        && ~isempty(strmatch('normal_y', fields)) ...
        && ~isempty(strmatch('normal_z', fields)))
    % Read the normal vectors.
    normal_x = data(:, strmatch('normal_x', fields));
    normal_y = data(:, strmatch('normal_y', fields));
    normal_z = data(:, strmatch('normal_z', fields));
    
    % Organize the values, if the point cloud is organized.
    if (height > 1)
        normal_x = reshape(normal_x, height, width);
        normal_y = reshape(normal_y, height, width);
        normal_z = reshape(normal_z, height, width);
        normal = cat(3, normal_x, normal_y, normal_z);
    else
        normal = [normal_x, normal_y, normal_z];
    end
end

%% Read DAT file.
% Create the DAT file's name.
datFilename = [filename(1 : end-length('.pcd')), '_info.dat'];

% Verify that the file exists.
fid = fopen(datFilename, 'r');
position = [];
if (fid ~= -1)
    position = posread(datFilename);
end
       
%% Create the return structure.
outcloud.viewpoint = viewpoint;
outcloud.intensity = intensity;

if (~isempty(position))
    outcloud.time = position.time;
    outcloud.odometry = position.odometry;
    outcloud.gps = position.gps;
end

if (isempty(color) && isempty(normal))
    outcloud.pointCloud = pointCloud(location);
elseif (~isempty(color) && isempty(normal))
    outcloud.pointCloud = pointCloud(location, 'Color', color);
elseif (isempty(color) && ~isempty(normal))
    outcloud.pointCloud = pointCloud(location, 'Normal', normal);
else
    outcloud.pointCloud = pointCloud(location, ...
                                     'Color', color, 'Normal', normal);
end

end
