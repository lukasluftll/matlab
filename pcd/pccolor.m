function pcout = pccolor(pcin, color)
% PCCOLOR Colorize point cloud.
%   PCOUT = PCCOLOR(PCIN) assigns to each point of the point cloud PCIN a
%   color that encodes its z-coordinate using the current colormap. 
%   The darkest color is assigned to the lowest point, the brightest color 
%   is assigned to the highest point.
%
%   PCIN and PCOUT are pointCloud objects.
%
%   PCOUT = PCCOLOR(PCIN, COLOR) uses the real-valued HEIGHTxWIDTH matrix 
%   COLOR to colorize the point cloud, where HEIGHT and WIDTH are the
%   height and width of the point cloud.
%   Its lowest value corresponds to the darkest color, its highest value 
%   corresponds to the brightest color.
%
%   Example:
%      pc = pcread('teapot.ply');
%      pcshow(pccolor(pc, pc.Location(:,1)))
%
%   See also POINTCLOUD, COLORMAP.

% Copyright 2016 Alexander Schaefer

%% Validate input and output.
% Check the number of output and input arguments.
nargoutchk(0, 1)
narginchk(1, 2)

% Check the input point cloud.
validateattributes(pcin, {'pointCloud'}, {'scalar'}, '', 'PCIN')

% If COLOR is not specified, set it to the z-coordinates of the points.
if nargin < 2
    if ismatrix(pcin.Location) == 2
        color = pcin.Location(:,3);
    else
        color = pcin.Location(:,:,3);
    end
end

% Check COLOR.
pcsize = size(pcin.Location);
validateattributes(color, {'numeric'}, ...
    {'real', 'size', [pcsize(1:end-1),1]}, '', 'COLOR')

%% Colorize.
% Convert the real color values to indices in [1; 64].
color = floor(normm(color) * 63.9999) + 1;

% Transform the color indices to RGB colors using the colormap.
fig = figure('visible', 'off');
palette = uint8(colormap .* 255);
close(fig);

% Assign the RGB colors to the output point cloud.
pcout = pcin;
pcout.Color = reshape(palette(color(:),:), size(pcout.Location));

end
