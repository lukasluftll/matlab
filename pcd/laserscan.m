classdef laserscan < handle
    % LASERSCAN Object for storing 3D laser scan.
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS, RLIM, TFORM) creates a 
    %   laser scan object.
    %
    %   AZIMUTH, ELEVATION, and RADIUS are HEIGHTxWIDTH matrices defining
    %   the rays the laser sensor measured in spherical coordinates w.r.t.
    %   the sensor frame.
    %   AZIMUTH and ELEVATION must be finite angles in [rad].
    %   RADIUS may contain real and infinite values. All RADIUS values 
    %   outside interval RLIM are interpreted as no-return rays.
    %
    %   RLIM is a 2-element vector that defines the minimum and maximum 
    %   radius the sensor is able to measure.
    %
    %   TFORM is a homogeneous 4x4 transformation matrix defining the 
    %   transformation from a global reference frame to the laser sensor 
    %   frame. It defaults to identity.
    %
    %   Example:
    %      pcd = pcdread('castle.pcd');
    %      l = laserscan(pcd.azimuth, pcd.elevation, pcd.radius, [0, 100])
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
            % Check if AZIMUTH is a 2D matrix.
            if ~ismatrix(azimuth)
                error('AZIMUTH must be a 2D matrix.')
            end
            
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
            % Check if ELEVATION is a 2D matrix.
            if ~ismatrix(elevation)
                error('ELEVATION must be a 2D matrix.')
            end
            
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
            % Check if RADIUS is a 2D matrix.
            if ~ismatrix(radius)
                error('RADIUS must be a 2D matrix.')
            end
            
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
            % Check if TFORM is a homogeneous transformation matrix.
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
            
            % Check the all values are finite.
            if ~all(isfinite(rlim))
                error('RLIM must be finite.')
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
        % Construct laserscan object.
        function obj = laserscan(azimuth, elevation, radius, rlim, tform)
            % Check number of input arguments.
            narginchk(4, 5);
            
            % If the sensor pose is not given, set it to identity.
            if nargin < 5
                tform = eye(4);
            end
            
            % Store the input.
            obj.azimuth = azimuth;
            obj.elevation = elevation;
            obj.radius = radius;
            obj.tform = tform;
            obj.rlim = rlim;
        end
        
        % Identify the return and no-return rays.
        function r = ret(obj)
            % RET Identify return and no-return rays.
            %   R = RET(OBJ) returns for each ray of the laser scan 
            %   whether it was reflected by the environment or not.
            %
            %   R is a HEIGHTxWIDTH logical matrix, where HEIGHT and WIDTH 
            %   define the size of the spherical coordinate matrices
            %   AZIMUTH, ELEVATION, and RADIUS used when constructing the
            %   laserscan object. All true elements of IR correspond to
            %   returned rays, all false elements correspond to no-return
            %   rays.
            
            r = obj.radius >= obj.rlim(1) & obj.radius <= obj.rlim(2);
        end
        
        % Plot the laser scan.
        function plot(obj)
            % PLOT Plot laser scan.
            %   PLOT(OBJ) visualizes the rays originating from the laser 
            %   scanner. Returned rays are plotted in red, no-return rays 
            %   are plotted in light gray.

            %% Convert spherical to Cartesian coordinates.
            % Define the plotted length of all rays.
            plotradius = obj.radius;
            plotradius(~ret(obj)) = obj.rlim(2);
            
            % Convert the spherical coordinates to Cartesian coordinates.
            [x, y, z] = sph2cart(obj.azimuth(:), obj.elevation(:), ...
                plotradius(:));
            
            % Rotate the rays around the origin.
            p = zeros(numel(x), 3);
            for i = 1 : numel(x)
                p(i,:) = (obj.tform(1:3,1:3) * [x(i); y(i); z(i)])';
            end
            
            % Extract the scan origin from the transformation matrix.
            o = obj.tform(1:3,4)';

            %% Plot rays.
            % Plot the returned rays.
            r = ret(obj);
            pr = kron(p(r(:),:), [0; 1]);
            pr = pr + repmat(o, size(pr, 1), 1);
            rray = plot3(pr(:,1), pr(:,2), pr(:,3), 'Color', 'red');
            if ~isempty(rray)
                % Set transparency.
                rray.Color(4) = 0.5;
            end

            % Plot the no-return rays.
            pnr = kron(p(~r(:),:), [0; 1]);
            pnr = pnr + repmat(o, size(pnr, 1), 1);
            hold on
            nrray = plot3(pnr(:,1), pnr(:,2), pnr(:,3), 'Color', 'k');
            if ~isempty(nrray)
                % Set transparency.
                nrray.Color(4) = 0.03;
            end
            
            %% Plot ray endpoints.
            pcshow(pointCloud(pr), 'MarkerSize', 25);

            %% Plot decoration.
            % Plot the sensor origin.
            plot3(o(1), o(2), o(3), 'Color', 'k', ...
                'Marker', '.', 'MarkerSize', 50);

            % Plot the Cartesian axes.
            plot3(o(1) + [0, max(x(:))], o(2) + [0, 0], o(3) + [0, 0], ...
                'Color', 'r', 'LineWidth', 3);
            plot3(o(1) + [0, 0], o(2) + [0, max(y(:))], o(3) + [0, 0], ...
                'Color', 'g', 'LineWidth', 3);
            plot3(o(1) + [0, 0], o(2) + [0, 0], o(3) + [0, max(z(:))], ...
                'Color', 'b', 'LineWidth', 3);

            % Label the axes.
            xlabel('x'); ylabel('y'); zlabel('z');

            axis equal
            grid on
            hold off
        end
    end
end
