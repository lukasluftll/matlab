function ray = pc2sph(point, model)
% PC2SPH Convert a point cloud to spherical coordinates.
%   RAY = PC2SPH(POINT) takes an MxNx3 matrix that contains an organized
%   point cloud in Cartesian coordinates and converts it to spherical 
%   coordinates. 
%   RAY is the resulting MxNx3 matrix.
%   RAY(:,:,1) are the azimuth angles in radians,
%   RAY(:,:,2) are the elevation angles in radians, 
%   RAY(:,:,3) are the radius values.
%
%   RAY = PC2SPH(POINT, MODEL) uses Lidar sensor-specific information to
%   reconstruct the azimuth and elevation angles of Cartesian NaN points.
%   The radius of a point whose Cartesian coordinates are NaN is set to 
%   Inf.
%   The available models are:
%
%      'generic'   - (default) same behavior as without MODEL argument
%      'vlp16'     - Velodyne VLP-16; can handle single and dual returns
%
%   See also CART2SPH, NAN, INF.

% Copyright 2016 Alexander Schaefer

%% Validate input.
% Check if the user provided enough input arguments.
narginchk(1, 2)

% Check if the input point cloud is organized.
if size(point, 3) ~= 3 || ndims(point) ~= 3
    error('Input point cloud must be of size MxNx3].')
end

% If no sensor model argument is provided, set it to default.
if nargin < 2
    model = 'generic';
end

%% Convert coordinates.
[a, e, r] = cart2sph(point(:,:,1), point(:,:,2), point(:,:,3));

%% Perform model-specific angle reconstruction.
switch model
    case 'vlp16'
        % Check if the rows are properly organized.
        esign = sign(diff(e));
        eorg = (all(esign(~isnan(esign))) == 1 ...
            || all(esign(~isnan(esign))) == -1);
        if ~eorg
            error('Point cloud must be properly organized.')
        end
        
        % Determine whether the point cloud contains single or dual
        % returns.
        i = repmat([1 0], size(a, 1), size(a, 2)/2);
        singleReturn = all(a(i) ~= a(~i));
        a = a(i);
        e = e(i);
        
        % Make sure the azimuth angles increase along the rows.
        for col = 2 : size(a, 2)
            a(:,col) = a(:,col) + 2*pi * (a(:,col)<a(:,col-1));
        end

        % Reconstruct missing elevation values.
        em = repmat(nanmean(e, 2), 1, size(e, 2));
        e(isnan(e)) = em(isnan(e));
        
        % Reconstruct missing azimuth values.
        for row = 1 : size(a, 1)
            ari = 1 : size(a, 2);
            ar = a(row,:);
            ar(isnan(ar)) = interp1(ari(~isnan(ar)), ar(~isnan(ar)), ...
                ari(isnan(ar)), 'linear', 'extrap');
            a(row,:) = ar;
        end
        
        % Crop azimuth angles to (-pi; pi].
        a(a > pi) = a(a > pi) - 2*pi;
        
        % Expand azimuth and elevation matrices for dual return 
        % point clouds.
        if singleReturn
            a = reshape([a; a], size(a, 1), []);
            e = reshape([e; e], size(e, 1), []);
        end
        
        % Set radius of NaN points to Inf.
        r(isnan(r)) = Inf;
end

%% Construct return matrix.
ray = cat(3, a, cat(3, e, r));

end
