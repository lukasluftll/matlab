function pcdwrite(pcd, file)
% PCDWRITE Write point cloud to PCD file.
%   PCDWRITE(PCD, FILE) writes the point cloud data PCD to the
%   ASCII-coded PCD file FILE.
%
%   PCD is either a pointCloud object or a struct with an arbitrary number 
%   of fields.
%
%   If PCD is a struct, the size of the payload data matrices contained in
%   the fields is HEIGHTxWIDTHxCOUNT. HEIGHT and WIDTH indicate the height 
%   and the width of the point cloud, whereas COUNT specifies the dimension 
%   of the data of one point. COUNT can vary from matrix to matrix.
%
%   Some field names are predefined and may not be used for payload data:
%      - VERSION,
%      - FIELDS,
%      - SIZE,
%      - TYPE,
%      - COUNT,
%      - WIDTH,
%      - HEIGHT,
%      - VIEWPOINT,
%      - POINTS,
%      - DATA.
%
%   Except for PCD.VIEWPOINT, all these fields are ignored. However,
%   PCD.VIEWPOINT indicates the pose of the sensor frame relative to the 
%   reference frame of all points. PCD.VIEWPOINT must be a 4x4 homogeneous 
%   transformation matrix. If PCD.VIEWPOINT is not defined, the VIEWPOINT 
%   header entry of the PCD file is set to identity.
%
%   FILE is the name of the PCD file. It may or may not contain the file
%   extension.
%
%   Example:
%      pcdwrite(pcread('teapot.ply'), 'teapot.pcd')
%
%   See also PCWRITE, PCDREAD.
 
% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check the number of input arguments.
narginchk(2, 2)

% Check the input arguments types.
validateattributes(pcd, {'struct'}, {'numel', 1}, '', 'PCD')
validateattributes(file, {'char'}, {'row', 'nonempty'}, '', 'FILE')

% If the file is given without extension, append the extension.
ext = '.pcd';
if numel(file) <= numel(ext) || ~strcmpi(file(end-numel(ext)+1:end), ext)
    file = [file, ext];
end

%% Preprocess input.
% Extract the viewpoint.
field = fieldnames(pcd);
vp = strcmpi('viewpoint', field);
if isempty(find(vp, 1))
    viewpoint = eye(4);
elseif ishrt(pcd.(field{vp}))
    viewpoint = pcd.(field{vp});
else
    error('PCD.VIEWPOINT must be a 4x4 homogeneous transformation.')
end

% Find all predefined fields in PCD.
rem = false(size(field));
prefield = {'version', 'fields', 'size', 'type', 'count', 'width', ...
    'height', 'viewpoint', 'points', 'data'};
for i = 1 : numel(prefield)
    rem = rem | strcmpi(prefield{i}, field);
end
field = field(~rem);

% Check the content of the remaining fields of PCD.
for i = 1 : numel(field)
    if ~isnumeric(pcd.(field{i}))
        error(['PCD.', toupper(field{i}), ...
            ' does not contain numeric values.'])
    end
    if ndims(pcd.(field{1})) ~= ndims(pcd.(field{i})) ...
            || any(size(pcd.(field{1})) ~= size(pcd.(field{i})))
        error(['Sizes of PCD.', toupper(field{1}), ' and ', ...
            toupper(field{i}), ' do not match.'])
    end
end

%% Open file.
% Try to open the file.
fid = fopen(file, 'w');
if fid == -1
    error(['Cannot create file ', file, '.'])
end

% Make sure the file is closed upon function termination.
cleaner = onCleanup(@() fclose(fid));

%% Write PCD header.
% Write the opening line.
fprintf(fid, '# .PCD v.7 - Point Cloud Data file format\n');

% PCD file version.
fprintf(fid, 'VERSION .7\n');

% Fields contained in the PCD file.
fprintf(fid, 'FIELDS ');
cellfun(@(f) fprintf(fid, '%s ', f), field);
fprintf(fid, '\b\n');

% Sizes of the datatypes.
fprintf(fid, 'SIZE ');
for i = 1 : numel(field)
    fielddata = pcd.(field{i});
    datainfo = whos('fielddata');
    fprintf(fid, '%i ', datainfo.bytes / numel(fielddata));
end
fprintf(fid, '\b\n');

% Names of the datatypes.
fprintf(fid, 'TYPE ');
    function t = pcdtype(d)
        % PCDTYPE PCD datatype code letter of given data.
        switch class(d)
            case {'single', 'double'}
                t = 'F';
            case {'int8', 'int16', 'int32', 'int64'}
                t = 'I';
            case {'uint8', 'uint16', 'uint32', 'uint64'}
                t = 'U';
            otherwise
                t = '?';
        end
    end
cellfun(@(f) fprintf(fid, '%s ', pcdtype(pcd.(f))), field);
fprintf(fid, '\b\n');

% Number of elements per field per point.
fprintf(fid, 'COUNT ');
cellfun(@(f) fprintf(fid, '%i ', size(pcd.(f), 3)), field);
fprintf(fid, '\b\n');

% Point cloud width.
width = size(pcd.(field{1}),2);
fprintf(fid, 'WIDTH %i\n', width);

% Point cloud height.
height = size(pcd.(field{1}),1);
fprintf(fid, 'HEIGHT %i\n', height);

% Acquisition viewpoint of the points in the dataset.
fprintf(fid, 'VIEWPOINT ');
vpv = [tform2trvec(viewpoint), tform2quat(viewpoint)];
arrayfun(@(c) fprintf(fid, '%i ', c), vpv);
fprintf(fid, '\b\n');

% Total number of points.
fprintf(fid, 'POINTS %i\n', width*height);

% Data type of the PCD file.
fprintf(fid, 'PCD ascii\n');

%% Write payload data.
% Compile the data matrix.
dlmdata = struct2cell(pcd).';
dlmdata = dlmdata(~rem);
dlmdata = cellfun(@(d) permute(d,[2,3,1]), dlmdata, 'UniformOutput',false);
dlmdata = cellfun(@(d) reshape(d,[],size(d,2),1), dlmdata, ...
    'UniformOutput',false);
dlmdata = cell2mat(dlmdata);

% Write the payload data to file.
dlmwrite(file, dlmdata, 'delimiter',' ', 'precision','%.7g', '-append');
fprintf('\n');

end
