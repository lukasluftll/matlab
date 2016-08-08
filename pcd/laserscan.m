classdef laserscan < handle
    % LASERSCAN Object for storing 3D laser scan.
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS, RLIM, TFORM) creates a 
    %   laser scan object.
    %
    %   AZIMUTH, ELEVATION, and RADIUS are N-element vectors defining
    %   the rays the laser sensor measured in spherical coordinates w.r.t.
    %   the sensor frame. N is the number of measured rays.
    %   AZIMUTH and ELEVATION must be finite angles in [rad].
    %   RADIUS may contain real and infinite values. All RADIUS values 
    %   outside interval RLIM are interpreted as no-return rays.
    %
    %   RLIM is a 2-element vector that defines the minimum and maximum 
    %   radius the sensor is able to measure. It defaults to the minimum
    %   and maximum given RADIUS.
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
            % Check if AZIMUTH contains the same number of elements as the 
            % current azimuth data.
            if ~isempty(obj.azimuth) && numel(azimuth)~=numel(obj.azimuth)
                error(['AZIMUTH must contain ', ...
                    num2str(numel(obj.azimuth)), 'elements.'])
            end
            
            % Make sure all angles are finite.
            if ~all(isfinite(azimuth(:)))
                error('AZIMUTH must be finite.')
            end
            
            % Assign new azimuth data to object.
            obj.azimuth = azimuth(:);
        end
        
        % Change the elevation angles.
        function set.elevation(obj, elevation)            
            % Check if ELEVATION contains the same number of elements as 
            % the current elevation data.
            if ~isempty(obj.elevation) ...
                    && numel(elevation)~=numel(obj.elevation)
                error(['ELEVATION must contain ', 
                    num2str(numel(obj.elevation)), 'elements.'])
            end
            
            % Make sure all angles are finite.
            if ~all(isfinite(elevation(:)))
                error('ELEVATION must be finite.')
            end
            
            % Assign new elevation data to object.
            obj.elevation = elevation(:);
        end
        
        % Change the ray radii.
        function set.radius(obj, radius)           
            % Check if RADIUS contains the size as the current radius data.
            if ~isempty(obj.radius) && numel(radius)~=numel(obj.radius)
                error(['RADIUS must contain ', ...
                    num2str(numel(obj.radius)), 'elements.'])
            end
            
            % Set all infinite values to NaN.
            radius(~isfinite(radius)) = NaN;
            
            % Assign new radius data to object.
            obj.radius = radius(:);
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
            narginchk(3, 5);
            
            % If the sensor range is not given, set it to the minimum and
            % maximum RADIUS.
            if nargin < 4
                rlim = [min(radius(:)), max(radius(:))];
            end
            
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
        
        % Return number of rays.
        function n = count(obj)
            % Return number of rays, both returned rays and no-returns.
            
            n = numel(obj.azimuth);
        end
        
        % Get Cartesian ray direction vectors.
        function v = cart(obj)
            % CART(OBJ) Get Cartesian ray direction vectors.
            %   V = CART(OBJ) returns a Nx3 matrix that contains the
            %   Cartesian direction vectors of the rays of the laser scan
            %   w.r.t. the global reference frame. N is the number of rays.
            %   The length of the direction vectors of no-return rays is
            %   set to unity.
            
            % Convert spherical to Cartesian coordinates.
            finiteradius = obj.radius;
            finiteradius(~ret(obj)) = 1;
            [x, y, z] = sph2cart(obj.azimuth, obj.elevation, finiteradius);
                        
            % Rotate the ray vectors around the origin.
            v = zeros(numel(x), 3);
            for i = 1 : numel(x)
                v(i,:) = (obj.rotation * [x(i); y(i); z(i)])';
            end
        end
        
        % Identify the return and no-return rays.
        function r = ret(obj)
            % RET Identify return and no-return rays.
            %   R = RET(OBJ) returns for each ray of the laser scan 
            %   whether it was reflected by the environment or not.
            %
            %   R is a N-element logical vector. All true elements of R 
            %   correspond to returned rays, all false elements correspond 
            %   to no-return rays.
            
            r = obj.radius >= obj.rlim(1) & obj.radius <= obj.rlim(2);
        end
        
        % Get the sensor position.
        function p = position(obj)
            % POSITION Laser sensor position.
            %   P = POSITION(OBJ) returns a 1x3 vector that specifies the 
            %   laser sensor position.
            
            p = obj.tform(1:3,4)';
        end
        
        % Get the sensor orientation.
        function r = rotation(obj)
            % ROTATION Laser sensor orientation.
            %   R = ROTATION(OBJ) returns a 3x3 rotation matrix that 
            %   specifies the laser sensor orientation.
            
            r = obj.tform(1:3,1:3);
        end
        
        % Plot the laser scan.
        function plot(obj)
            % PLOT Plot laser scan.
            %   PLOT(OBJ) visualizes the rays originating from the laser 
            %   scanner. Returned rays are plotted in red, no-return rays 
            %   are plotted in light gray.

            %% Prepare data.
            % Get the logical indices of the returned rays.
            ir = ret(obj);
            
            % Convert spherical to Cartesian coordinates.
            p = cart(obj);
            
            % Set the plotted length of the no-return rays to the maximum
            % sensor range.
            p(~ir,:) = p(~ir,:) * obj.rlim(2);
            
            % Get the sensor position.
            o = obj.position;

            %% Plot rays.
            % Plot the returned rays.
            pr = kron(p(ir(:),:), [0; 1]);
            pr = pr + repmat(o, size(pr, 1), 1);
            lsr = plot3(pr(:,1), pr(:,2), pr(:,3), 'Color', 'red');
            if ~isempty(lsr)
                % Set transparency.
                lsr.Color(4) = 0.5;
            end

            % Plot the no-return rays.
            pnr = kron(p(~ir(:),:), [0; 1]);
            pnr = pnr + repmat(o, size(pnr, 1), 1);
            hold on
            lsnr = plot3(pnr(:,1), pnr(:,2), pnr(:,3), 'Color', 'k');
            if ~isempty(lsnr)
                % Set transparency.
                lsnr.Color(4) = 0.03;
            end
            
            %% Plot ray endpoints.
            pcshow(pointCloud(pr), 'MarkerSize', 25);

            %% Plot decoration.
            % Plot the sensor origin.
            plot3(o(1), o(2), o(3), 'Color', 'k', ...
                'Marker', '.', 'MarkerSize', 50);
            
            % Plot the Cartesian axes.
            plot3(o(1) + [0,max(p(:,1))], o(2) + [0,0], o(3) + [0,0], ...
                'Color', 'r', 'LineWidth', 3);
            plot3(o(1) + [0,0], o(2) + [0,max(p(:,2))], o(3) + [0,0], ...
                'Color', 'g', 'LineWidth', 3);
            plot3(o(1) + [0,0], o(2) + [0,0], o(3) + [0,max(p(:,3))], ...
                'Color', 'b', 'LineWidth', 3);

            % Label the axes.
            xlabel('x'); ylabel('y'); zlabel('z');

            axis equal
            grid on
            hold off
        end
    end
end
