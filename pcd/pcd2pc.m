function pc = pcd2pc(pcd)

nargoutchk(1, 1)
narginchk(1, 1)

validateattributes(pcd, {'struct'}, {}, '', 'PCD')

if ~isfield(pcd, 'x') || ~isfield(pcd, 'y') || ~isfield(pcd, 'z')
    error('PCD does not contain fields X, Y, and Z.')
end

sizechk(pcd.x, pcd.y, pcd.z);
location = cat(3, pcd.x, cat(3, pcd.y, pcd.z));
pc = pointCloud(location);

field = fieldnames(pcd);
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
