function ray = pc2sph(point)
% PC2SPH Convert a point cloud to spherical coordinates.
%   RAY = PC2SPH(POINT) takes an MxNx3 matrix that contains an 
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
%   elevation angles. The range coordinate of the point is set to Inf.
%
%   See also CART2SPH, NAN, INF.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided enough input arguments.
narginchk(1, 1)

% Check if the input point cloud is organized.
if size(point, 3) ~= 3 || ndims(point) ~= 3
    error('Input point cloud must be of size MxNx3].')
end

%% Convert coordinates.
% Convert the Cartesian coordinates to spherical coordinates.
[a, e, r] = cart2sph(point(:,:,1), point(:,:,2), point(:,:,3));

%% Reconstruct missing azimuth values.
% Determine the logical indices of the point columns whose azimuth values 
% cross the -pi to +pi border.
altcol = repmat(any(sign(a)==1) & any(sign(a)==-1), size(a, 1), 1);

% Make sure the azimuth values within one column do not cross the -pi to 
% +pi border. This would compromise the computation of the mean.
a(altcol & a<0) = a(altcol & a<0) + 2*pi;

% Calculate the mean azimuth values over the columns.
acm = nanmean(a);

% Reconstruct the missing values.
i = 1 : length(acm);
acm(isnan(acm)) = interp1(i(~isnan(acm)), acm(~isnan(acm)), ... 
    i(isnan(acm)), 'linear', 'extrap');
acm = repmat(acm, size(a, 1), 1);

% Calculate the mean offset of the azimuth values from the mean over 
% the rows.
arm = repmat(nanmean(a-acm, 2), 1, size(a, 2));

% Create an approximated azimuth matrix using the mean and the mean offset.
am = acm + arm;
am(am > pi) = am(am > pi) - 2*pi;

% Reconstruct missing azimuth values.
a(isnan(a)) = am(isnan(a));

%% Reconstruct missing elevation angles.
% Create an elevation angle matrix from the elevation mean values.
erm = nanmean(e, 2);

% Interpolate or extrapolate missing values.
erm(isnan(erm)) = interp1((1:length(erm)) .* ~isnan(erm), ...
    erm(~isnan(erm)), erm(isnan(erm)), 'linear', 'extrap');
em = repmat(erm, 1, size(e, 2));

% Reconstruct missing elevation angles from the matrix created by
% interpolation and extrapolation.
e(isnan(e)) = em(isnan(e));

%% Define missing radius values.
% NaN points are reflected at infinity, so set the radius to Inf.
r(isnan(r)) = Inf;

%% Construct return matrix.
ray = cat(3, a, cat(3, e, r));

end
