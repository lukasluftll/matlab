classdef laserscan
    % LASERSCAN Object for storing 3D laser scan.
    %   L = LASERSCAN(SP, AZIMUTH, ELEVATION, RADIUS, RLIM) creates a 
    %   laser scan object.
    %
    %   SP is a 4x4xN matrix, where N is the number of rays of the laser 
    %   scan. Its n-th page defines the laser sensor pose at the 
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
    %   radius that the sensor is able to measure. It defaults to the 
    %   minimum and maximum given RADIUS.
    %
    %   LASERSCAN properties:
    %   SP         - Sensor poses in global reference frame
    %   AZIMUTH    - Azimuth angles in sensor frame [rad]
    %   ELEVATION  - Elevation angles in sensor frame [rad]
    %   RADIUS     - Ray lengths
    %   RLIM       - Sensor measurement range
    %
    %   LASERSCAN methods:
    %   COUNT      - Number of rays
    %   END2CART   - Cartesian coordinates of ray endpoints in global frame
    %   DIR2CART   - Cartesian ray direction vectors in global frame
    %   LS2PC      - Transform laser scan to point cloud
    %   RET        - Identify returned rays
    %   LSCONCAT   - Concatenate laser scans
    %   CHREF      - Change global reference system
    %   SELECT     - Select subset of measurements
    %   PLOT       - Plot laser scan
    %
    %   See also TRAV, SLAB.
    
    % Copyright 2016 Alexander Schaefer
    
    properties
        % SP 4x4xN matrix; concatenation of 4x4 homogeneous transformation 
        % matrices defining for each ray the transformation from the 
        % global reference frame of the laser scan to the laser sensor.
        % N is the number of rays.
        sp;
        
        % AZIMUTH N-element vector defining the azimuth angle of the rays 
        % with respect to the sensor frame in [rad]. The azimuth angle is 
        % the counterclockwise angle in the x-y plane measured from the 
        % positive x axis. N is the number of measured rays.
        azimuth;
        
        % ELEVATION N-element vector defining the elevation angle of the  
        % rays with respect to the sensor frame in [rad]. The elevation 
        % angle is the angle measured from the x-y plane. N is the number 
        % of measured rays.
        elevation;
        
        % RADIUS N-element vector defining the length of the rays
        % originating from the laser sensor. N is the number of measured 
        % rays.
        % NaN values identify no-return rays.
        radius;
        
        % RLIM 2-element vector defining the minimum and maximum radius the
        % sensor is able to measure.
        rlim;
    end
    
    methods
        function obj = set.sp(obj, sp)
            % SET.SP Set sensor poses.
            
            % Check if SP contains the same number of dimensions and
            % elements as the current sensor pose data.
            if ~isempty(obj.sp) && (ndims(sp) ~= ndims(obj.sp) ...
                    || ~all(size(sp) == size(obj.sp)))
                error(['SP must be a ', num2str(size(obj.sp,1)), 'x', ...
                    num2str(size(obj.sp,2)), 'x', ...
                    num2str(size(obj.sp,3)), ' matrix.'])
            end
            
            % Check if SP contains homogeneous transformation matrices 
            % only.
            if ~all(ishrt(sp))
                error(['SP(:,:,', num2str(find(~ishrt(sp), 1)), ...
                        ') is not a homogeneous transformation.'])
            end
            
            % Assign new sensor pose data to object.
            obj.sp = sp;
        end
        
        function obj = set.azimuth(obj, azimuth)
            % SET.AZIMUTH Set azimuth angles in sensor frame in [rad].
            
            % Check if AZIMUTH contains the same number of elements as the 
            % current azimuth data.
            if ~isempty(obj.azimuth) && numel(azimuth)~=numel(obj.azimuth)
                error(['AZIMUTH must contain ', ...
                    num2str(numel(obj.azimuth)), ' elements.'])
            end
            
            % Make sure all angles are finite.
            if ~all(isfinite(azimuth(:)))
                error('AZIMUTH must be finite.')
            end
            
            % Assign new azimuth data to object.
            obj.azimuth = azimuth(:);
        end
        
        function obj = set.elevation(obj, elevation)
            % SET.ELEVATION Set elevation angles in sensor frame in [rad].
            
            % Check if ELEVATION contains the same number of elements as 
            % the current elevation data.
            if ~isempty(obj.elevation) ...
                    && numel(elevation)~=numel(obj.elevation)
                error(['ELEVATION must contain ', ...
                    num2str(numel(obj.elevation)), ' elements.'])
            end
            
            % Make sure all angles are finite.
            if ~all(isfinite(elevation(:)))
                error('ELEVATION must be finite.')
            end
            
            % Assign new elevation data to object.
            obj.elevation = elevation(:);
        end
        
        function obj = set.radius(obj, radius)
            % SET.RADIUS Set the ray radii.
            
            % Check if RADIUS contains the size as the current radius data.
            if ~isempty(obj.radius) && numel(radius)~=numel(obj.radius)
                error(['RADIUS must contain ', ...
                    num2str(numel(obj.radius)), ' elements.'])
            end
            
            % Set all infinite values to NaN.
            radius(~isfinite(radius)) = NaN;
            
            % Assign new radius data to object.
            obj.radius = radius(:);
        end
        
        function obj = set.rlim(obj, rlim)
            % SET.RLIM Set sensor radius limits.
            
            %% Validate input.
            % Check the number of input arguments.
            narginchk(2, 2)
            
            % Check the size of the limit vector.
            if numel(rlim) ~= 2
                error('RLIM must contain exactly 2 elements.')
            end
            
            % Check if the minimum limit is positive.
            if rlim(1) < 0
                error('RLIM(1) must be greater or equal zero.')
            end
            
            % Check if all values are finite.
            if ~all(isfinite(rlim))
                error('RLIM must be finite.')
            end
            
            % Check if the limit vector is ordered.
            if diff(rlim) <= 0
                error('RLIM(2) must be greater than RLIM(1).');
            end
            
            %% Assign new limits.
            obj.rlim = rlim;
        end
        
        function idxchk(obj, i)
            % IDXCHK Check laser scan ray indices.
            %   IDXCHK(OBJ, I) throws an error if I is not a valid index
            %   matrix for the rays of the laserscan object OBJ.
            
            if islogical(i)
                if numel(i) > obj.count
                    error('I contains too many elements.')
                end
            elseif isnumeric(i)
                if any(rem(i, 1) ~= 0)
                    error('I must contain positive interger values only.')
                end
                if any(i < 1 | i > obj.count)
                    error('I contains invalid indices.')
                end
            else
                error('I must be a logical or numeric vector.')
            end
        end
    end
    
    methods ( Access = public )
        function obj = laserscan(sp, azimuth, elevation, radius, rlim)
            % LASERSCAN Constructor.
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(4, 5);
            
            % If the sensor range is not given, set it to the minimum and
            % maximum RADIUS.
            if nargin < 5
                radiusFinite = radius(isfinite(radius));
                rlim = [min(radiusFinite), max(radiusFinite)];
            end
            
            % Check whether the sizes of the input arguments match.
            n = size(sp,3);
            if numel(azimuth) ~= n || numel(elevation) ~= n ...
                    || numel(radius) ~= n
                error('Number of elements of input arguments must match.')
            end
            
            %% Assign input.
            obj.sp = sp;
            obj.azimuth = azimuth;
            obj.elevation = elevation;
            obj.radius = radius;
            obj.rlim = rlim;
        end
        
        function n = count(obj)
            % COUNT Number of rays.
            n = numel(obj.azimuth);
        end
        
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
            
            % If no index is given, select all rays.
            if nargin < 2
                i = 1 : obj.count;
            end
            
            % Check validity of indices.
            obj.idxchk(i)
            
            % If the laser scan or the index vector is empty, abort.
            p = zeros(0, 3);
            if obj.count < 1 || isempty(i) || all(~i)    
                return
            end
            
            %% Compute coordinates.
            % Compute the ray endpoint coordinates in the sensor frame.
            [x,y,z] = sph2cart(obj.azimuth(i(:)), obj.elevation(i(:)), ...
                obj.radius(i(:)));
           
            % Transform the ray endpoints into the reference frame of the
            % laser scan.
            p = reshape(cart2hom([x,y,z]).', 4, 1, []);
            p = hom2cart(permute(pagetimes(obj.sp(:,:,i(:)), p), [3,1,2]));
        end
        
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
            
            % Check validity of indices.
            obj.idxchk(i)
            
            % If the laser scan or the index vector is empty, return.
            v = zeros(0, 3);
            if obj.count < 1 || isempty(i) || all(~i)    
                return
            end
            
            %% Compute direction vectors.
            % Compute the normalized ray direction vectors in the sensor 
            % frame.
            [x,y,z] = sph2cart(obj.azimuth(i(:)), obj.elevation(i(:)), 1);
            
            % Transform the ray direction vectors into the reference frame
            % of the laser scan.
            rot = rotm2tform(tform2rotm(obj.sp(:,:,i(:))));
            v = hom2cart(pagetimes(rot, cart2hom([x,y,z]).').');
        end
        
        function pc = ls2pc(obj)
            % LS2PC Laser scan to point cloud.
            %   PC = LS2PC(OBJ) returns a pointCloud object that contains 
            %   the ray endpoints of the laser scan.
            %   
            %   No-return rays of the laser scan result in NaN points.
            %
            %   The reference frames of the point cloud and of the laser 
            %   scan are the same.   
            pc = pointCloud(end2cart(obj));
        end
        
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
            
            % Create a matrix whose rows contain the sensor ranges of the
            % individual scans.
            xrlim = reshape([x.rlim], 2, []);
            
            % Build new laserscan object.
            y = laserscan(ysp, ...
                reshape([x.azimuth], [], 1), ...
                reshape([x.elevation], [], 1), ...
                reshape([x.radius], [], 1), ...
                [max(xrlim(1,:)), min(xrlim(2,:))]);
        end
        
        function obj = chref(obj, tform)
            % CHREF Change reference system.
            %   LS = CHREF(OBJ, TFORM) changes the reference system of the
            %   laserscan object OBJ according to the transformation TFORM.
            %   After this operation, the sensor poses SP are specified in
            %   the new reference frame.
            %
            %   TFORM is the 4x4 homogeous transformation matrix from the 
            %   current reference frame of the scan to the new reference 
            %   frame.
            
            %% Validate input.
            % Check the number of input arguments.
            narginchk(1, 2)
            
            % Check the validity of the homogeneous transformation matrix.
            hrtchk(tform)
            
            %% Change sensor poses.
            obj.sp=pagetimes(repmat(inv(tform),1,1,size(obj.sp,3)),obj.sp);
        end
        
        function sub = select(obj, i)
            % SELECT Select subset of measurements.
            %   SUB = SELECT(OBJ, I) returns a laserscan object that
            %   contains only the laser measurement indexed by I.
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
            
            %% Select subset of measurements.
            sub = laserscan(obj.sp(:,:,i), ...
                obj.azimuth(i), obj.elevation(i), obj.radius(i), obj.rlim);
        end
       
        function plot(obj)
            % PLOT Plot laser scan.
            %   PLOT(OBJ) visualizes the rays originating from the laser 
            %   scanner. Returned rays are plotted in red. No-return rays 
            %   are plotted in light gray with length equal to maximum 
            %   sensor range.
            %
            %   If the scan contains a large set of rays, a random 
            %   selection of rays and sensor poses is plotted to minimize
            %   computational effort.
          
            %% Select subset of rays to plot.           
            % Build logical index vector that selects the maximum number of
            % rays at random.
            i = false(obj.count, 1);
            i(randsample(1:obj.count, min([obj.count, 2000]))) = true;
            
            %% Compute sensor positions and ray endpoints.
            % Identify returned and no-return rays to plot.
            ir = ret(obj) & i;
            inr = ~ret(obj) & i;
            
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
            
            % Get the sensor origin for each ray.
            sr = ht2tv(lsr.sp);
            snr = ht2tv(lsnr.sp);
            
            %% Plot.
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
            
            % Plot sensor poses.
            l = 0.05 * mean(sqrt(sum(tform2trvec(obj.sp).^2, 2)));
            i = round(linspace(1, obj.count, min([10, obj.count])));
            plotht(obj.sp(:,:,i), l);
            ptext = tform2trvec(obj.sp(:,:,1));
            text(ptext(1), ptext(2), ptext(3), 'Sensor pose')
            
            % Plot decoration.
            plotht(eye(4), min(max(er)), 'LineWidth', 2)
            text(0, 0, 0, 'Scan reference')
            labelaxes
            grid on
        end
    end
end
