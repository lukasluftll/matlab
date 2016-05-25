function outcloud = pccolor(incloud, palette, interval)
% pccolor Colorize a 3D point cloud using intensity values.
%   outcloud = pccolor(incloud, palette, interval) uses a 64-step color
%   palette to transform the given intensity values into colors and 
%   writes them to the point cloud object. 
%   
%   The input parameters must have the following properties:
%
%   incloud     structure that contains at least the following elements:
%               pointCloud  - pointCloud object, organized or unorganized.
%               intensity   - matrix that specifies intensity values.
%                             For unorganized point clouds, 
%                             its size must be [height * width, 1]; 
%                             for organized point clouds, its size must be
%                             [height, width, 1].
%
%   palette     optional MATLAB colormap name. Frequently used colormaps
%               are 'gray', 'hot', and 'parula'. Type 'help colormap' 
%               for more palettes.
%
%               Default: 'jet'.
%
%   interval    2-element vector that contains the intensity values that
%               correspond to the darkest and the brightest color.
%               Intensity values below interval(1) or above interval(2) 
%               will be set to interval(1) and interval(2), respectively.
%
%               Default: [0; 100].
%
%   The return value is described as follows:
%
%   outcloud    pointCloud object with color values that correspond to the 
%               given intensities.
%
%   Example : colorize a point cloud
%   --------------------------------
%   cloud = pcdread('terrain.pcd');
%   cloud = pccolor(cloud, 'hot', [0; 255]);
%   pcshow(cloud);
%
%   See also pointCloud, pcshow, colormap, pcdread.
 
% Copyright 2016 Alexander Schaefer

%% Check input.
% If the user does not specify the required arguments, abort.
if (nargin < 1)
    error('Invalid input: not enough input arguments.')
end

% If the user does not specify a color map, use the default color map.
if (nargin < 2)
    palette = 'jet';
end

% If the user does not specify an interval, use the default interval.
if (nargin < 3)
    interval = [0; 100];

% Get the dimensions of the point cloud.
pcsize = size(incloud.pointCloud.Location);

% Check whether the size of the intensity matrix matches the size of the
% point cloud.
intensity = incloud.intensity;
if (size(pcsize) >= 3)
    if (size(intensity) ~= pcsize(1:2))
        error(['Invalid input: ', ...
            'intensity matrix size does not match point cloud size.']);
    end
else
    if (size(intensity, 1) ~= pcsize(1))
        error(['Invalid input: ', ...
            'intensity matrix size does not match point cloud size.']);
    end
end

% Check whether the given interval vector is valid.
interval = interval(1:2);
if (interval(1) > interval(2))
    warning('Invalid input: interval(1) > interval(2).');
end
if (~isreal(interval))
    warning('Invalid input: interval is not real.');
end

% Clip the intensities to the given interval.
intensity(intensity < interval(1)) = interval(1);
intensity(intensity > interval(2)) = interval(2);

% Scale the intensities to [0; 63].
intensity = round((intensity - interval(1)) * 63 / diff(interval));

% Invert the intensities, as high intensities correspond to bright colors,
% and fit them to [1; 64].
intensity = 64 - intensity;

%% Colorize point cloud.
% Transform the colormap, which contains 64 RGB values in [0; 1.0], to
% RGB values in [0; 255], as required for coloring point clouds.
% Function colormap requires an open figure.
colormapFigure = figure; 
palette = uint8(colormap(palette) .* 255);
close(colormapFigure);

% Convert the intensities to colors.
color = palette(intensity(:), :);

% Reshape the color matrix, so it is organized if the point cloud 
% is organized.
color = reshape(color, pcsize);

% Copy the manipulated color matrix to the output cloud.
outcloud = copy(incloud.pointCloud);
outcloud.Color = color;

end
