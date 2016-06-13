function outcloud = ptsread(filename, subsampling)
% PTSREAD Read a 3D point cloud from PTS file.
%   OUTCLOUD = PCDREAD(FILENAME) reads a point cloud from the 
%   PTS file specified by the string filename.
%   OUTCLOUD = PCDREAD(FILENAME, SUBSAMPLING) reads every 
%   SUBSAMPLING-th line of the given PTS file.
%
%   The return value OUTCLOUD is a structure with the following elements:
%
%   pointCloud      pointCloud object (unorganized point cloud)
%
%   intensity       intensity values of the point cloud as column vector
%
%   PTS file contents
%   -----------------
%   PTS files contain the following information: 
%   * 3D coordinates of the points
%   * scalar intensity of each point
%   * RGB color value of each point
%
%   Example : read a point cloud from a PTS file
%   --------------------------------------------
%   outcloud = ptsread('field.pts');
%   pcshow(outcloud.pointCloud);
%
%   See also PCREAD, PCDREAD, POINTCLOUD, PCWRITE, PCSHOW.
 
%  Copyright 2016 Alexander Schaefer

%% Validate the input.
% Make sure the user specified the correct number of arguments.
narginchk(1, 2);

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

% If there is no subsampling factor given, set it to 1.
if nargin < 2
    subsampling = 1;
end

%% Read the data.
% Read the number of points the file is assumed to contain.
count = textscan(fgetl(fid), '%d');
count = count{1};
fclose(fid);

% Read the data.
try
    data = dlmread(filename, ' ', 1, 0);
catch exception
    warning([exception.message, ' Stopped reading PTS file.']);
end

% Check the file's validity.
if count ~= size(data, 1)
    warning('PTS file header specifies incorrect number of points.');
end

% Read the point locations.
location = data(1 : subsampling : size(data, 1), 1 : 3);

% Read the intensities.
outcloud.intensity = [];
if size(data, 2) >= 4
    outcloud.intensity = data(1 : subsampling : size(data, 1), 4);
end

% Read the color values.
if size(data, 2) >= 7
    color = uint8(data(1 : subsampling : size(data, 1), 5 : 7));
end
       
%% Create the point cloud.
if isempty(color)
    outcloud.pointCloud = pointCloud(location);
else
    outcloud.pointCloud = pointCloud(location, 'Color', color);
end

end
