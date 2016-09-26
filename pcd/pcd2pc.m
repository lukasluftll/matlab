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

nargoutchk(0, 1)
narginchk(1, 1)

validateattributes(pcd, {'struct'}, {}, '', 'PCD')

field = fieldnames(pcd);
ix = find(strcmpi('x', field), 1);
iy = find(strcmpi('y', field), 1);
iz = find(strcmpi('z', field), 1);
if isempty(ix) || isempty(iy) || isempty(iz)
    error('PCD does not contain fields X, Y, and Z.')
end

x = pcd.(field{ix});
y = pcd.(field{iy});
z = pcd.(field{iz});
sizechk(x, y, z);
location = cat(3, x, cat(3, y, z));
pc = pointCloud(location);


ic = find(strcmpi('rgb', field) | strcmpi('color', field), 1);
if ~isempty(ic)
    color = pcd.(field{ic});
    sizechk(location, color)
    pc.Color = color;
end

ii = find(strcmpi('intensity', field) | strcmpi('intensities', field), 1);
if isempty(ic) && ~isempty(ii)
    intensity = pcd.(field{ii});
    sizechk(location(:,:,1), intensity)
    
    intensity = floor((intensity-min(intensity(:))) * 63.9999 / range(intensity)) + 1;
    
    % Transform the colormap, which contains 64 RGB values in [0; 1.0], to
    % RGB values in [0; 255], as required for coloring point clouds.
    % Function colormap requires an open figure.
    fig = figure('visible', 'off'); 
    palette = uint8(colormap('jet') .* 255);
    close(fig);

    % Convert the intensities to colors.
    pc.Color = reshape(palette(intensity(:),:), size(location));
end

end
