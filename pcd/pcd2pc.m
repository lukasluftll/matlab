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
location = cat(3, x, cat(3, y, z));
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
        % Check the size of the intensity matrix.
        intensity = pcd.(field{ii});
        sizechk(location(:,:,1), intensity)
    
        % Convert the intensities into indices in [1; 64].
        intensity = (intensity-min(intensity(:))) / range(intensity);
        intensity = floor(normm(intensity) * 63.9999) + 1;
    
        % Transform the intensity values to colors using the colormap.
        fig = figure('visible', 'off'); 
        palette = uint8(colormap('jet') .* 255);
        close(fig);
        pc.Color = reshape(palette(intensity(:),:), size(location));
    end 
else
    % Check the size of the color matrix.
    color = pcd.(field{ic});
    sizechk(location, color)
    
    % Add color information to the pointCloud object.
    pc.Color = color;
end

end
