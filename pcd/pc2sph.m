function beams = pc2sph(points)
% pc2sph Convert a point cloud to spherical coordinates.
%   beams = pc2sph(points) takes an [M x N x 3] matrix that contains an 
%   organized point cloud and converts the points to spherical coordinates. 
%   The columns of the returned [M x N x 3] matrix beams contain the 
%   azimuth, elevation, and range of the points.
%
%   Conversion of NaN points
%   ------------------------
%   If the input matrix contains NaN points, pc2sph uses the point 
%   coordinates around each NaN point to interpolate its azimuth and 
%   elevation angles. Its range is set to NaN.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided enough input arguments.
if nargin < 1
    error('Not enough input arguments.');
end

% Check if the input point cloud is organized.
if size(points, 3) ~= 3
    error('Input point cloud must be of size [M x N x 3].');
end

%% Convert coordinates.
% Convert the Cartesian coordinates to spherical coordinates.
p = reshape(points, numel(points(:,:,1)), 3);
[azimuth, elevation, radius] = cart2sph(p(:,1), p(:,2), p(:,3));
beams = reshape([azimuth, elevation, radius], size(points));

% Compute the mean azimuth angles of the point columns and the mean 
% elevation angles of the point rows.
azimuth = nanmean(beams(:,:,1), 1);
elevation = nanmean(beams(:,:,2), 2);

% If a row or a colums contains NaN values only, interpolate the 
% corresponding azimuth or elevation angle using the mean angles 
% computed above.
azimuth(isnan(azimuth)) ...
    = interp1(1:length(azimuth), azimuth, isnan(azimuth));
elevation(isnan(elevation)) ...
    = interp1(1:length(elevation), elevation, isnan(elevation));

% Set the azimuth and elevation angles of all NaN points.
for row = 1 : size(points, 1)
    for col = 1 : size(points, 2)
        if isnan(beams(row,col,:))
            beams(row,col,1:2) = [azimuth(col), elevation(row)];
        end
    end
end

end
