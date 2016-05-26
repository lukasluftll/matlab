function groundIndex = coneseg(incloud, coneAngle)
% CONESEG Segments the ground in the given point cloud.
%   GROUNDINDEX = CONESEG(INCLOUD) puts a cone that opens along +z 
%   on top of each point. All points that are not contained in the cone 
%   of any other point are segmented as ground points.
%
%   GROUNDINDEX = CONESEG(INCLOUD, CONEANGLE) does not use the default 
%   cone opening angle of sqrt(2) rad, but the value specified by the user.
%
%   Example : segment and display ground points 
%   -------------------------------------------
%   cloud = pcdread('terrain.pcd');
%   cloud = cloud.pointCloud;
%   pcshow(cloud.select(coneseg(cloud)));
%
%   See also POINTCLOUD, PCSHOW.
 
%  Copyright 2016 Alexander Schaefer

%% Prepare the input data.
% If the user did not specify a cone opening angle, do it now.
if (nargin < 2)
    coneAngle = pi/4;
end

% Check the validity of the cone opening angle input argument.
coneAngle = abs(coneAngle);
if (coneAngle > pi/2)
    coneAngle = pi/2;
    warning('Invalid input: cone opening angle > pi/2 set to pi/2.');
end
cosConeAngle = cos(coneAngle);

% If the input cloud is organized, reshape the point location matrix from
% 3D to 2D for easier computation.
if (size(incloud.Location, 3) > 1)
    location = [reshape(incloud.Location(:, :, 1), [incloud.Count, 1]), ...
        reshape(incloud.Location(:, :, 2), [incloud.Count, 1]), ...
        reshape(incloud.Location(:, :, 3), [incloud.Count, 1])];
else
    location = incloud.Location;
end

% Store the maximum z coordinate of all points.
zMax = max(incloud.ZLimits);

% Create a vector that stores a value for each point of the cloud.
% The value greater than 0 indicates that the point has been  
% segmented as ground.
isGround = 1 : size(location, 1);

% Create a progress bar.
waitbarHandle = waitbar(0, 'Segmenting ground ...');

%% Perform the cone test.
% Iterate through all points and perform the cone test.
progress = 0.0;
for (p = 1 : size(location, 1)) %#ok<*NO4LP>
    % If the point has already been detected to fall into a cone, skip it.
    if (isGround(p) > 0)    
        % Store the coordinates of the current point.
        point = location(p, :);

        % Calculate the size of the cone.
        coneHeight = zMax - point(3);
        coneRadius = coneHeight * tan(coneAngle);

        % Calculate the region of interest for this point.
        roi = [point(1) - coneRadius, point(1) + coneRadius; 
            point(2) - coneRadius, point(2) + coneRadius;
            point(3), zMax];

        % Get the points in the region of interest.
        roiIndex = findPointsInROI(incloud, roi);
        reshape(roiIndex, [numel(roiIndex), 1]);
        roiPoints = location(roiIndex, :);

        % Calculate the vectors originating from the cone tip to the
        % points of the region of interest.
        d = roiPoints - repmat(point, size(roiPoints, 1), 1);

        % Check whether or not the points fall into the cone.
        cosAlpha = d(:, 3) ./ sqrt(sum(d.^2, 2));
        nongroundIndex = roiIndex(cosAlpha >= cosConeAngle);
        isGround(nongroundIndex) = 0;
    end

    % Advance the process bar every 1/100.
    if (p / size(location, 1) > progress + 0.01)
        progress = p / size(location, 1);
        waitbar(progress, waitbarHandle);
    end
end

%% Post-process the result.
% If the input point cloud is organized, reshape the segmentation result
% matrix from 1D to 2D.
if (size(incloud.Location, 3) > 1)
    isGround = reshape(isGround, ...
        [size(incloud.Location, 1), size(incloud.Location, 2)]);
end

% Remove the indices of the points segmented as non-ground.
isGround(isGround <= 0) = [];
groundIndex = isGround;

% Close the process bar.
close(waitbarHandle);

end

