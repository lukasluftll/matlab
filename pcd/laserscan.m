classdef laserscan < handle
    % LASERSCAN Object for storing 3D laser scan.
    %   L = LASERSCAN(SP, AZIMUTH, ELEVATION, RADIUS, RLIM) creates a 
    %   laser scan object.
    %
    %   SP is a 4x4xN matrix, where N is the number of rays the laser 
    %   scan contains. Its n-th page defines the laser sensor pose at the 
    %   moment of the n-th observation; it is a 4x4 homogeneous
    %   transformation matrix from the reference frame of the scan to the 
    %   laser sensor pose.
    %
    %   AZIMUTH, ELEVATION, and RADIUS are N-element vectors defining
    %   the rays the laser sensor measured in spherical coordinates with
    %   respect to the sensor frame.
    %   AZIMUTH and ELEVATION must be finite angles in [rad].
    %   RADIUS may contain real and infinite values. All RADIUS values 
    %   outside interval RLIM are interpreted as no-return rays.
    %
    %   RLIM is a 2-element vector that defines the minimum and maximum 
    %   radius the sensor is able to measure. It defaults to the minimum
    %   and maximum given RADIUS.
    %
    %   LASERSCAN properties:
    %   SP         - Sensor poses in global reference frame
    %   AZIMUTH    - Azimuth angles [rad]
    %   ELEVATION  - Elevation angles [rad]
    %   RADIUS     - Ray lengths
    %   RLIM       - Sensor measurement range
    %
    %   LASERSCAN methods:
    %   COUNT      - Number of rays
    %   END2CART   - Cartesian coordinates of ray endpoints
    %   DIR2CART   - Cartesian ray direction vectors
    %   LS2PC      - Transform laser scan to point cloud
    %   RET        - Identify returned rays
    %   LSCONCAT   - Concatenate laser scans
    %   PLOT       - Plot laser scan
    %
    %   See also TRAV, SLAB.
    
    % Copyright 2016 Alexander Schaefer
    
    properties
        % SP 4x4xN matrix; concatenation of 4x4 homogeneous transformation 
        % matrices defining for each ray the transformation from the 
        % reference frame of the laser scan to the laser sensor.
        % N is the number of rays.
        sp;
        
        % AZIMUTH N-element vector defining the azimuth angle of the rays 
        % with respect to the sensor frame in [rad]. 
        % N is the number of measured rays.
        azimuth;
        
        % ELEVATION N-element vector defining the elevation angle of the  
        % rays with respect to the sensor frame in [rad]. 
        % N is the number of measured rays.
        elevation;
        
        % RADIUS N-element vector defining the length of the rays
        % originating from the laser sensor. 
        % N is the number of measured rays.
        % NaN values identify no-return rays.
        radius;
        
        % RLIM 2-element vector defining the minimum and maximum radius the
        % sensor is able to measure.
        rlim;
    end
    
    methods
        % Change the sensor poses.
        function set.sp(obj, sp)
            % Check if SP contains the same number of dimensions and
            % elements as the current sensor pose data.
            if ~isempty(obj.sp)
                if ndims(sp) ~= ndims(obj.sp)
                    error(['SP must be a ', num2str(ndims(obj.sp)), ...
                        'D matrix.'])
                end
                if any(size(sp) ~= size(obj.sp))
                    error(['SP must be of size ', ...
                        num2str(size(obj.sp, 1)), 'x', ...
                        num2str(size(obj.sp, 2)), 'x', ...
                        num2str(size(obj.sp, 3)), '.'])
                end
            end
            
            % Check if SP contains homogeneous transformation matrices.
            hrt = ishrt(sp);
            if ~all(hrt)
                error(['SP(:,:,', num2str(find(hrt, 1)), ...
                        ') is not a homogeneous transformation.'])
            end
            
            % Assign new sensor pose data to object.
            obj.sp = sp;
        end
        
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
                error(['ELEVATION must contain ', ...
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
        function obj = laserscan(sp, azimuth, elevation, radius, rlim)
            % Check number of input arguments.
            narginchk(4, 5);
            
            % If the sensor range is not given, set it to the minimum and
            % maximum RADIUS.
            if nargin < 5
                radiusFinite = radius(isfinite(radius));
                rlim = [min(radiusFinite), max(radiusFinite)];
            end
            
            % Store the input.
            obj.sp = sp;
            obj.azimuth = azimuth;
            obj.elevation = elevation;
            obj.radius = radius;
            obj.rlim = rlim;
        end
        
        % Return number of rays.
        function n = count(obj)
            % COUNT Number of rays.
            n = numel(obj.azimuth);
        end
        
        % Get Cartesian coordinates of ray endpoints.
        function p = end2cart(obj, i)
            % END2CART(OBJ) Cartesian coordinates of ray endpoints.
            %   P = LS2CART(OBJ) returns an Nx3 matrix that contains the
            %   Cartesian coordinates of the ray endpoints with respect to
            %   the reference frame of the laser scan.
            %   N is the number of rays.
            %
            %   P = LS2CART(OBJ, I) returns the Cartesian coordinates of
            %   the endpoints of selected rays. I is an M-element vector
            %   that contains the indices of the rays of interest.
            %   P is a Mx3 matrix that contains the Cartesian endpoints of
            %   these rays.
            %
            %   The coordinates of no-return rays are set to NaN.
            
            %% Validate input.
            narginchk(1, 2)
            
            % If no rays are selected, select all.
            if nargin < 2
                i = 1 : obj.count;
            end
            
            %% Compute coordinates.
            % Compute the ray endpoint coordinates in the sensor frame.
            [x,y,z] = sph2cart(obj.azimuth(i(:)), obj.elevation(i(:)), ...
                obj.radius(i(:)));
           
            % Transform the ray endpoints into the reference frame of the
            % laser scan.
            p = hom2cart(httimes(obj.sp(:,:,i(:)), cart2hom([x,y,z]).').');
        end
        
        % Get Cartesian ray direction vectors.
        function v = dir2cart(obj, i)
            % DIR2CART(OBJ) Cartesian ray direction vectors.
            %   V = DIRECTION2CART(OBJ) returns an Nx3 matrix that contains
            %   the normalized Cartesian direction vectors of the rays with
            %   respect to the reference frame of the laser scan.
            %   N is the number of rays.
            %
            %   DIR2CART(OBJ, I) returns the normalized Cartesian direction 
            %   vectors of selected rays. I is an M-element vector that
            %   contains the indices of the rays of interest.
            %   P is a Mx3 matrix that contains the direction vectors of
            %   these rays.
            
            %% Validate input.
            narginchk(1, 2)
            
            % If no rays are selected, select all.
            if nargin < 2
                i = 1 : obj.count;
            end
            
            %% Compute direction vectors.
            % Compute the ray direction vectors in the sensor frame.
            [x,y,z] = sph2cart(obj.azimuth(i(:)), obj.elevation(i(:)), 1);
            
            % Transform the ray direction vectors into the reference frame
            % of the laser scan.
            rot = rotm2tform(tform2rotm(obj.sp(:,:,i(:))));
            v = hom2cart(httimes(rot, cart2hom([x,y,z]).').');
        end
        
        % Transform laser scan to point cloud.
        function pc = ls2pc(obj)
            % LS2PC Transform laser scan to point cloud.
            %   PC = LS2PC(OBJ) returns a pointCloud object that contains 
            %   the ray endpoints of the laser scan.
            %   
            %   No-return rays of the laser scan result in NaN points.
            %
            %   The reference frames of the point cloud and of the laser 
            %   scan are the same.   
            pc = pointCloud(end2cart(obj));
        end
        
        % Identify the return and no-return rays.
        function r = ret(obj)
            % RET Identify return and no-return rays.
            %   R = RET(OBJ) returns for each ray of the laser scan 
            %   whether it was reflected by the environment or not.
            %
            %   R is an N-element logical vector. All true elements of R 
            %   correspond to returned rays, all false elements correspond 
            %   to no-return rays.
            r = obj.radius >= obj.rlim(1) & obj.radius <= obj.rlim(2);
        end
        
        % Concatenate laser scans.
        function y = lsconcat(x)
            % LSCONCAT Concatenate laser scans.
            %   Y = LSCONCAT(X) concatenates all laserscan objects 
            %   contained in the laserscan vector X to build a laserscan
            %   object Y.
            %
            %   The sensor range of Y is the common subset of the ranges of
            %   all laserscan objects in X.
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(1, 1)
            
            % Check type of input argument.
            if ~isa(x, 'laserscan')
                error('X must be a laserscan object.')
            end
            
            %% Concatenate data.
            % Concatenate sensor poses.
            ysp = x(1).sp; 
            for i = 2 : numel(x)
                ysp = cat(3, ysp, x(i).sp);
            end
            
            % Create matrix that contains sensor ranges.
            xrlim = reshape([x.rlim], 2, []);
            
            % Build new laserscan object.
            y = laserscan(ysp, reshape([x.azimuth], [], 1), ...
                reshape([x.elevation], [], 1), ...
                reshape([x.radius], [], 1), ...
                [max(xrlim(1,:)), min(xrlim(2,:))]);
        end            
       
        % Plot the laser scan.
        function plot(obj)
            % PLOT Plot laser scan.
            %   PLOT(OBJ) visualizes the rays originating from the laser 
            %   scanner. Returned rays are plotted in red, no-return rays 
            %   are plotted in light gray. 
          
            % Limit number of rays to plot.
            nl = 2000;
            il = false(obj.count, 1);
            il(round(linspace(1, obj.count, min([obj.count, nl])))) = true;
            
            % Identify returned and no-return rays.
            ir = ret(obj) & il;
            inr = ~ret(obj) & il;
            
            % Create laserscan objects that contain the returned rays to
            % plot and the no-returns to plot.
            lsr = laserscan(obj.sp(:,:,ir), obj.azimuth(ir), ...
                obj.elevation(ir), obj.radius(ir), obj.rlim);            
            lsnr = laserscan(obj.sp(:,:,inr), obj.azimuth(inr), ...
                obj.elevation(inr), repmat(obj.rlim(2), sum(inr), 1), ...
                obj.rlim);
            
            % Compute plotted ray endpoints.
            er = end2cart(lsr);
            enr = end2cart(lsnr);
            
            % Get sensor origin for each ray.
            sr = tform2trvec(lsr.sp);
            snr = tform2trvec(lsnr.sp);
            
            % Plot returned rays.
            retplot = plot3([sr(:,1).'; er(:,1).'], ...
                [sr(:,2).'; er(:,2).'], ...
                [sr(:,3).'; er(:,3).']);
            set(retplot, 'Color', [1,0,0,0.3])
            
            % Plot no-return rays.
            hold on
            nrplot = plot3([snr(:,1).'; enr(:,1).'], ...
                [snr(:,2).'; enr(:,2).'], ...
                [snr(:,3).'; enr(:,3).']);
            set(nrplot, 'Color', [1,1,1,0.03]);
            
            % Plot point cloud.
            pcshow(pointCloud(er), 'MarkerSize', 80);
            
            % Plot sensor position.
            s = tform2trvec(obj.sp);
            plot3(s(:,1), s(:,2), s(:,3), 'k.');
            
            % Plot decoration.
            plotht(eye(4), min(max(er)), 'LineWidth', 2)
            labelaxes
            grid on
        end
    end
end
