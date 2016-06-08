function ray = pc2sph(points)
% PC2SPH Convert a point cloud to spherical coordinates.
%   RAY = PC2SPH(POINTS) takes an MxNx3 matrix that contains an 
%   organized point cloud in Cartesian coordinates and converts it 
%   to spherical coordinates. 
%   RAY is the resulting MxNx3 matrix.
%   RAY(:,:,1) are the azimuth angles in radians,
%   RAY(:,:,2) are the elevation angles in radians, 
%   RAY(:,:,3) are the range values.
%
%   Conversion of NaN points
%   ------------------------
%   If the input matrix contains NaN points, PC2SPH uses the point 
%   coordinates around each NaN point to interpolate its azimuth and 
%   elevation angles. Its range is set to NaN.
%
%   See also CART2SPH.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided the correct number of input arguments.
narginchk(1, 1);

% Check if the input point cloud is organized.
if size(points, 3) ~= 3 || ndims(points) ~= 3
    error('Input point cloud must be of size MxNx3].');
end

%% Convert coordinates.
% Convert the Cartesian coordinates to spherical coordinates.
p = reshape(points, numel(points(:,:,1)), 3);
[azimuth, elevation, radius] = cart2sph(p(:,1), p(:,2), p(:,3));
ray = reshape([azimuth, elevation, radius], size(points));

% Compute the mean azimuth angles of the point columns and the mean 
% elevation angles of the point rows.
azimuth = nanmean(ray(:,:,1), 1);
elevation = nanmean(ray(:,:,2), 2);

% If a row or a colums contains NaN values only, interpolate the 
% corresponding azimuth or elevation angle using the mean angles 
% computed above.
azimuthIndex = 1 : length(azimuth);
azimuth(isnan(azimuth)) ...
    = interp1(azimuthIndex, azimuth, azimuthIndex(isnan(azimuth)));
elevationIndex = 1 : length(elevation);
elevation(isnan(elevation)) ...
    = interp1(elevationIndex, elevation, elevationIndex(isnan(elevation)));

% Set the azimuth and elevation angles of all NaN points.
for row = 1 : size(ray, 1)
    for col = 1 : size(ray, 2)
        if isnan(ray(row,col,:))
            ray(row,col,1:2) = [azimuth(col), elevation(row)];
        end
    end
end

end
