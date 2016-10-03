function pc = pcd2pc(pcd)
% PCD2PC Convert struct to pointCloud object.
%   PC = PCD2PC(PCD) converts the point cloud data contained in the fields
%   of struct PCD to a pointCloud object. 
%
%   The recognized fields of PCD are:
%   X                      (required)  -  Cartesian x-coordinates
%   Y                      (required)  -  Cartesian y-coordinates
%   Z                      (required)  -  Cartesian z-coordinates
%   RGB/COLOR              (optional)  -  color
%   INTENSITY/INTENSITIES  (optional)  -  remission intensity
%
%   The X, Y, Z and INTENISTY matrices must be of size HEIGHTxWIDTH, where 
%   HEIGHT and WIDTH are the height and width of the point cloud.
%
%   RGB must be of size HEIGHTxWIDTHx3, where the pages contain the red,
%   green, and blue values respectively.
%
%   PC is a pointCloud object. If PCD specifies colors, this information
%   is copied to PC.COLOR. If PCD specifies no colors, but intensities, 
%   PC.COLOR contains the intensities translated into colors.
%
%   Example:
%      pcd = pcdread('pcd/data/office.pcd');
%      pcd2pc(pcd)
%
%   See also POINTCLOUD, PCDREAD.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check the number of input and output arguments.
nargoutchk(0, 1)
narginchk(1, 1)

% Check whether PCD is a struct.
validateattributes(pcd, {'struct'}, {}, '', 'PCD')

% Check whether PCD contains all required fields.
field = fieldnames(pcd);
ix = find(strcmpi('x', field), 1);
iy = find(strcmpi('y', field), 1);
iz = find(strcmpi('z', field), 1);
if isempty(ix) || isempty(iy) || isempty(iz)
    error('PCD does not contain fields X, Y, and Z.')
end

% Check the size of the Cartesian coordinate matrices.
validateattributes(pcd.(field{ix}), {'numeric'}, {'real'}, '', 'PCD.X')
validateattributes(pcd.(field{iy}), {'numeric'}, ...
    {'real', 'size', size(pcd.(field{ix}))}, '', 'PCD.Y')
validateattributes(pcd.(field{iz}), {'numeric'}, ...
    {'real', 'size', size(pcd.(field{ix}))}, '', 'PCD.Z')

%% Create pointCloud object.
location = cat(3,pcd.(field{ix}),cat(3,pcd.(field{iy}),pcd.(field{iz})));
pc = pointCloud(location);

%% Add color information.
% Look for color information in the struct.
ic = find(strcmpi('rgb', field) | strcmpi('color', field), 1);

% If the struct does not contain color information, use the intensities to
% color the point cloud.
if isempty(ic)
    % Check whether the struct contains intensity information.
    ii = find(strcmpi('intensity',field) | strcmpi('intensities',field),1);

    % If the struct contains intensity information, convert it to colors.
    if ~isempty(ii)
        % Check the intensity matrix.
        validateattributes(pcd.(field{ii}), {'numeric'}, ...
            {'real', 'size', size(location(:,:,1))}, ...
            '', ['PCD.', upper(field{ii})])
        
        % Convert the intensities to colors.
        pc = pccolor(pc, pcd.(field{ii}));
    end 
else
    % Check the color matrix.
    validateattributes(pcd.(field{ic}), {'uint8'}, ...
        {'size', size(pc.Location)}, '', ['PCD.', upper(field{ic})])
    
    % Add color information to the pointCloud object.
    pc.Color = pcd.(field{ic});
end

end
