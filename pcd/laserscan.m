classdef laserscan < handle
    % LASERSCAN Object for storing a 3D laser scan.
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS) creates a laser scan
    %   object.
    %
    %   AZIMUTH, ELEVATION, and RADIUS are HEIGHTxWIDTH matrices which
    %   define the rays the laser sensor measured in spherical coordinates.
    %   AZIMUTH and ELEVATION must be finite angles in [rad].
    %   RADIUS may contain infinite values, which are interpreted as 
    %   no-return rays.
    %
    %   The laser sensor pose is assumed to be equal to the reference 
    %   frame.
    %
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS, TFORM) additionally
    %   specifies the laser sensor pose with respect to the global 
    %   reference frame.
    %
    %   TFORM is a homogeneous 4x4 transformation matrix defining rotation 
    %   and translation only, no sheering, scaling, etc.
    %
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS, TFORM, RLIM) additionally
    %   defines the minimum and maximum radius the sensor is able to 
    %   measure using the ordered 2-element vector RLIM.
    %
    %   Example:
    %      pcd = pcdread('castle.pcd');
    %      l = laserscan(pcd.azimuth, pcd.elevation, pcd.radius)
    %
    %   See also TRAV, SLAB.
    
    % Copyright 2016 Alexander Schaefer
    
    properties
        tform;
        azimuth;
        elevation;
        radius;
        rlim;
    end
    
    methods
        % Change the azimuth angles.
        function set.azimuth(obj, azimuth)
            % Check if AZIMUTH has the same size as the current azimuth 
            % data.
            if ~isempty(obj.azimuth) && size(azimuth) ~= size(obj.azimuth)
                error(['AZIMUTH must be a matrix of size ', ...
                    num2str(size(obj.azimuth, 1)), 'x', ...
                    num2str(size(obj.azimuth, 2)), '.'])
            end
            
            % Make sure all angles are finite.
            if ~all(isfinite(azimuth(:)))
                error('AZIMUTH must be finite.')
            end
            
            % Assign new azimuth data to object.
            obj.azimuth = azimuth;
        end
        
        % Change the elevation angles.
        function set.elevation(obj, elevation)
            % Check if ELEVATION has the same size as the current elevation
            % data.
            if ~isempty(obj.elevation) ...
                    && size(elevation)~=size(obj.elevation)
                error(['ELEVATION must be a matrix of size ', ...
                    num2str(size(obj.elevation, 1)), 'x', ...
                    num2str(size(obj.elevation, 2)), '.'])
            end
            
            % Make sure all angles are finite.
            if ~all(isfinite(elevation(:)))
                error('ELEVATION must be finite.')
            end
            
            % Assign new elevation data to object.
            obj.elevation = elevation;
        end
        
        % Change the ray radii.
        function set.radius(obj, radius)
            % Check if RADIUS has the same size as the current radius data.
            if ~isempty(obj.radius) && size(radius) ~= size(obj.radius)
                error(['RADIUS must be a matrix of size ', ...
                    num2str(size(obj.radius, 1)), 'x', ...
                    num2str(size(obj.radius, 2)), '.'])
            end
            
            % Set all infinite values to NaN.
            radius(~isfinite(radius)) = NaN;
            
            % Assign new radius data to object.
            obj.radius = radius;
        end
        
        % Change the sensor pose.
        function set.tform(obj, tform)
            if ~ishrt(tform)
                error(['TFORM must be a homogeneous ', ...
                    '4x4 rotation-translation matrix.'])
            end
            
            % Assign new pose to object.
            obj.tform = tform;
        end
        
        % Change the sensor radius limits.
        function set.rlim(obj, rlim)
            % Check the size of the limit vector.
            if numel(rlim) ~= 2
                error('RLIM must contain exactly 2 elements.')
            end
            
            % Check the minimum limit is positive.
            if rlim(1) < 0
                error('RLIM(1) must be greater or equal zero.')
            end
            
            % Check the limit vector is ordered.
            if diff(rlim) <= 0
                error('RLIM(2) must be greater than RLIM(1).');
            end
            
            % Assign new limits to object.
            obj.rlim = rlim;
        end
    end
    
    methods ( Access = public )
        % Construct a laserscan object.
        function obj = laserscan(azimuth, elevation, radius, tform, rlim)
            % Check number of input arguments.
            narginchk(3, 5);

            % If the sensor pose is not given, set it to identity.
            if nargin < 4
                tform = eye(4);
            end
            
            % If the limit vector is not given, define it.
            if nargin < 5
                rlim = [0; +Inf];
            end
            
            % Store the input.
            obj.azimuth = azimuth;
            obj.elevation = elevation;
            obj.radius = radius;
            obj.tform = tform;
            obj.rlim = rlim;
        end
        
        % Identify the no-return rays.
        function inr = noret(obj)
            % NORET Identify no-return rays.
            %   INR = NORET(OBJ) returns the subscript indices of the
            %   no-return rays of the laser scan.
            %
            %   INR is a HEIGHTxWIDTH logical matrix, where HEIGHT and
            %   WIDTH define the size of the spherical coordinate matrices
            %   AZIMUTH, ELEVATION, and RADIUS used when constructing the
            %   laserscan object. All elements of INR that correspond to
            %   no-return rays are true, all others are false.
            
            inr = ~isfinite(obj.radius) | obj.radius < obj.rlim(1) ...
                | obj.radius > obj.rlim(2);
        end
        
        % Plot the laser scan.
        function plot(obj)
            % PLOT Plot laser scan.
            %   PLOT(OBJ) visualizes the rays originating from the laser 
            %   scanner. Returned rays are plotted in red, no-return rays 
            %   are plotted in light gray.
                        
            %% Compute visualized radius of no-return rays.
            plotradius = obj.radius;
            if isfinite(obj.rlim(2))
                plotradius(noret(obj)) = obj.rlim(2);
            else
                plotradius(noret(obj)) = max(obj.radius);    
            end

            % Convert the spherical coordinates to Cartesian coordinates.
            [x, y, z] = sph2cart(obj.azimuth, obj.elevation, plotradius);
            
            % Rotate the rays around the origin.
            
            
            p = [x, y, z];

%% Plot rays.
% Plot the returned rays.
pret = kron(p(ret(:),:), [0; 1]);
retray = plot3(pret(:,1), pret(:,2), pret(:,3), 'Color', 'red');
if ~isempty(retray)
    retray.Color(4) = 0.5;
end

% Plot the no-return rays.
pnan = kron(p(~ret(:),:), [0; 1]);
hold on
nanray = plot3(pnan(:,1), pnan(:,2), pnan(:,3), 'Color', 'k');
if ~isempty(nanray)
    nanray.Color(4) = 0.03;
end

%% Plot decoration.
% Plot the sensor origin.
plot3(0, 0, 0, 'Color', 'k', 'Marker', '.', 'MarkerSize', 50);

% Plot the Cartesian axes.
plot3([0, max(x(:))], [0, 0], [0, 0], 'Color', 'r', 'LineWidth', 3);
plot3([0, 0], [0, max(y(:))], [0, 0], 'Color', 'g', 'LineWidth', 3);
plot3([0, 0], [0, 0], [0, max(z(:))], 'Color', 'b', 'LineWidth', 3);

% Plot the point cloud.
pcshow(pointCloud(pret), 'MarkerSize', 30);
hold off

% Label the axes.
xlabel('x')
ylabel('y')
zlabel('z')

axis equal
grid on

end

        end
    end   
end
