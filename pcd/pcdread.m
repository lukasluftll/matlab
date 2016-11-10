function data = pcdread(file)
% PCDREAD Read contents of point cloud from PCD file.
%   DATA = PCDREAD(FILE) reads header and data fields of the PCD file FILE 
%   and returns them in struct DATA. 
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
%      PCD file name   <file>.pcd
%      DAT file name   <file>_info.dat
%  
%   Example:
%      data = pcdread('office.pcd');
%      cloud = pointCloud([data.x; data.y; data.z]');
%      pcshow(cloud);
%
%   See also POSREAD, PCREAD, PTSREAD.
 
% Copyright 2016 Alexander Schaefer
%
% This function implements the PCD file format as specified by the 
% Point Cloud Library:
% http://pointclouds.org/documentation/tutorials/pcd_file_format.php.

%% Validate input and output.
% Check the numbers of input and output arguments.
nargoutchk(0, 1)
narginchk(1, 1)

% Check the input argument.
validateattributes(file, {'char'}, {'row', 'nonempty'}, '', 'FILE')

% Verify the file exists.
if exist(file, 'file') ~= 2
    error(['File ', file, ' does not exist.'])
end

%% Open file.
% Try to open the file.
fid = fopen(file, 'r');
if fid == -1
    error(['Cannot read ', file, '.'])
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

%% Read PCD payload data.
num_fieldnames = numel(fieldcount);
switch pcdtype{1}
        case 'ascii'
            format = '';
            for j=1:num_fieldnames
                % Get the right format out of the fieldtype
                % and build the right format structure for textscan
                out = fieldtype(j);
                switch fieldtype(j)
                    case 'I', typ = 'd';
                    case 'U', typ = 'u';
                    case 'F', typ = 'f';
                end
                format = [format '%' typ num2str(fieldsize(j)*8)];
            end
            c = textscan(fid, format, count);
            c = c.';
            % Transform the cell array into a structure array.
            for index=1:length(c)
                c{index}=c{index}.';
            end
            % If the data type of a raw data column is not double,
            % convert it to the data type given in the PCD file header.
            uint = regexp(fieldtype', '[uU]');
            c(uint) = cellfun(@uint32, c(uint), 'UniformOutput',...
                false);
            int = regexp(fieldtype', '[iI]');
            c(int) = cellfun(@int32, c(int), 'UniformOutput', false);
            int = regexp(fieldtype', '[fF]');
            c(int) = cellfun(@double, c(int), 'UniformOutput', false);

            % prepare the cell data for return
            data = cell2struct(c, fieldname, 1); 

        case 'binary'
                startPos_fid = ftell(fid);
                points = zeros(3, count);
                for j=1:num_fieldnames
                    % map IUF -> int, uint, float
                    switch fieldtype(j)
                        case 'I'
                            fmt = 'int';
                        case 'U'
                            fmt = 'uint';
                        case 'F'
                            fmt = 'float';
                    end
                    format = ['*' fmt num2str(fieldsize(j)*8)];
                    % repositioning of the file indicator
                    fseek(fid, startPos_fid + sum(fieldsize(1:j-1)),...
                        'bof');
                    % file read with right position and format
                    data = fread(fid, [1 count], format,...
                        sum(fieldsize)-fieldsize(j));
                    % save each row in points in every iteration
                    points(j,:) = data;
                end
                % transform the point matrices in the cell structure  
                % needed for the given structure 
                c = cell(3,1);
                for index=1:num_fieldnames
                    c{index} = points(index, :);
                end
                % If the data type of a raw data column is not double, 
                % convert it to the data type given in the PCD file
                % header.
                uint = regexp(fieldtype', '[uU]');
                c(uint) = cellfun(@uint32, c(uint),...
                    'UniformOutput', false);
                int = regexp(fieldtype', '[iI]');
                c(int) = cellfun(@int32, c(int),...
                    'UniformOutput', false);
                int = regexp(fieldtype', '[fF]');
                c(int) = cellfun(@double, c(int),...
                    'UniformOutput', false);
                data = cell2struct(c, fieldname, 1); 
        otherwise
                error('unknown or not supported DATA mode: %s', pcdtype{1});
%% Check header field COUNT.
if count ~= size(data.(fieldname{1}), 1)
    warning('PCD header field COUNT is inconsistent with data.')
end
  

%% Add header data to return structure.
data.version = pcdversion;
data.width = width;
data.height = height;
data.viewpoint = eye(4);
data.viewpoint(1:3,:) = [quat2rotm(viewpoint(4:7)'), viewpoint(1:3)];
data.count = count;

%% Read DAT file.
% Determine the name of the DAT file.
datFile = [file(1:end-length('.pcd')), '_info.dat'];

% Check if the DAT file exists.
if exist(datFile, 'file') ~= 2
    return
end

% Check if the DAT file can be opened for read access.
datFid = fopen(datFile, 'r');
if datFid == -1
    return
end
fclose(datFid);

% Read position information if available.
pos = posread(datFile);

% Append position information to return structure.
posfield = fieldnames(pos);
for i = 1 : numel(posfield)
    data.(posfield{i}) = pos.(posfield{i});
end
end