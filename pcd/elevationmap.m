classdef elevationmap
    % ELEVATIONMAP Object for storing an elevation data grid.
    %   EM = ELEVATIONMAP(DATA, RES) creates an elevation map in the x-y 
    %   plane with grid resolution RES.
    %
    %   RES is a positive scalar.
    %
    %   If DATA is a pointCloud object, the elevation of each grid tile is
    %   set to the maximum z coordinate of all points of point cloud DATA 
    %   that belong to that tile.
    %
    %   If DATA is a MxN matrix, the matrix values define the elevation of
    %   the respective map tile.
    %
    %   Example:
    %      em = elevationmap(pcread('teapot.ply'), 0.1)
    %
    %   ELEVATIONMAP properties:
    %   SUPPORT    - Minimum coordinates of the map
    %   EXTENSION  - Extension of the map in numbers of map tiles
    %   RESOLUTION - Edge length of each tile
    %   ELEVATION  - Elevation data
    %
    %   ELEVATIONMAP methods:
    %   GETELEV    - Elevation at coordinate
    %   SETELEV    - Assign elevation at coordinate
    %   LIMITS     - Spatial extension of map
    %   MASK       - Compute mask of tiles corresponding to given points
    %   DIFF       - Height difference between 3D points and elevation map
    %   MATCH      - Compute distance between point cloud and elevation map
    %   FILLNAN    - Estimate elevation values of NaN tiles
    %   SETPRIOR   - Set NaN tiles to value
    %   PLOT       - Visualize elevation map
    %
    %   See also POINTCLOUD.
    
    % Copyright 2016 Alexander Schaefer
    
    properties ( SetAccess = private )
        % SUPPORT Minimum coordinates of the map.
        %   [xmin, ymin].
        support
        
        % EXTENSION Extension of the map in numbers of map tiles.
        %   [nx, ny], with nd being the number of tiles in d-direction.
        extension
        
        % RESOLUTION Edge length of each tile.
        %   Scalar.
        resolution
    end
    
    properties
        % ELEVATION Elevation data.
        %   IxJ matrix with I and J being the number of tiles in x- and 
        %   y-direction.
        elevation
    end
    
    methods
        function obj = set.elevation(obj, elevation)
            % SET.ELEVATION Set elevation values.
            
            %% Validate input.
            narginchk(2, 2)
            if isempty(obj.elevation)
                validateattributes(elevation, {'numeric'}, ...
                    {'real', 'ndims', 2}, '', 'ELEVATION')
            else
                validateattributes(elevation, {'numeric'}, ...
                    {'real', 'size', size(obj.elevation)}, '', 'ELEVATION')
            end
            
            %% Assign values.
            obj.elevation = elevation;
        end
    end
    
    methods ( Access = private )
        function i = idx(obj, p)
            % IDX Compute tile index given x-y coordinates.
            %   I = IDX(OBJ, P) computes the linear index I of the tile to
            %   which the given x-y coordinates P correspond.
            %
            %   P is an Mx2 vector.
            %
            %   I is an M-element column vector. If P(m,:) lies outside the
            %   map, I(m) is NaN.
            
            %% Validate input and output.
            nargoutchk(0, 1)
            narginchk(2, 2)

            %% Compute tile index.
            % Determine which coordinates lie inside the map.
            lim = obj.limits;
            valid = lim(1,1) <= p(:,1) & p(:,1) < lim(1,2) ...
                & lim(2,1) <= p(:,2) & p(:,2) < lim(2,2);
            
            % Initialize the index vector.
            i = NaN(size(p,1), 1);
            
            % Compute the indices corresponding to the valid coordinates.
            % Make sure the indices stay within the valid range: Due to 
            % limited numerical precision, FLOOR() of a quotient may return
            % too large a value.
            d = p(valid,:) - repmat(obj.support, sum(valid), 1);
            isub = floor(d / obj.resolution) + 1;
            ix = constrain(isub(:,1), [1,obj.extension(1)]);
            iy = constrain(isub(:,2), [1,obj.extension(2)]);
            i(valid) = sub2ind(obj.extension, ix, iy);
        end                
    end
    
    methods ( Access = public )      
        function obj = elevationmap(data, res)
            % ELEVATIONMAP Constructor.
            
            %% Validate input and output.
            % Check number of input and output arguments.
            nargoutchk(0, 1)
            narginchk(2, 2)
            
            % Validate data.
            validateattributes(data, {'pointCloud', 'numeric'}, {'2d'}, ...
                '', 'DATA')
            
            % Check the validity of the point cloud.
            if isa(data, 'pointCloud')
                validateattributes(data, {'pointCloud'}, {'scalar'}, ...
                    '', 'DATA')
                validateattributes(data.XLimits, {'numeric'}, ...
                    {'finite', 'numel', 2}, '', 'PC.XLIMITS')
                validateattributes(data.YLimits, {'numeric'}, ...
                    {'finite', 'numel', 2}, '', 'PC.YLIMITS')
            end

            % Check the validity of the resolution.
            validateattributes(res, {'numeric'}, ...
                {'scalar', 'nonnegative'}, '', 'RES')

            %% Compute property values.
            % Store the resolution.
            obj.resolution = res;
            
            % Compute the property values if given the data directly.
            if isnumeric(data)
                obj.elevation = data;
                obj.support = [0,0];
                obj.extension = size(data);
                return
            end
            
            % Compute the property values if given a point cloud.
            % Compute the minimum x and y coordinates of the map if given
            % a point cloud.
            obj.support = floor([data.XLimits(1),data.YLimits(1)]/res)*res;
            
            % Compute the number of map tiles in x and y direction.
            obj.extension = floor(...
                ([data.XLimits(2),data.YLimits(2)]-obj.support) / res) + 1;
            
            % Remove all NaN and infinite points from the point cloud.
            data = removeInvalidPoints(data);
            
            % Make sure the point cloud is represented by an Mx3 vector.
            point = reshape(data.Location(:), data.Count, 3, 1);
            
            % Compute the indices of the tiles to which each point belongs.
            ie = obj.idx(point(:,1:2));
            
            % Set the elevation of each tile to the maximum z coordinate of
            % all points belonging to this tile.
            e = NaN(obj.extension);
            for i = 1 : size(point,1)
                e(ie(i)) = max(e(ie(i)), point(i,3));
            end
            obj.elevation = e;
        end
        
        function e = getelev(obj, c)
            % GETELEV Elevation at coordinate.
            %   E = GETELEV(OBJ, C) returns the elevation value 
            %   at coordinate C.
            %
            %   C is an Mx2 matrix. Each row specifies an x-y coordinate.
            %
            %   E is an M-element column vector. The m-th row contains the
            %   elevation value of the tile at coordinate C(m,:).
            %
            %   If C(m,:) lies outside the map, E(m) yields NaN.
            
            %% Validate input and output arguments.
            nargoutchk(0, 1)
            narginchk(2, 2)
            validateattributes(c, {'numeric'}, {'real','ncols',2}, '', 'C')
            
            %% Determine elevation at coordinate.
            i = obj.idx(c);
            e = NaN(size(c,1), 1);
            e(isfinite(i)) = obj.elevation(i(isfinite(i)));
        end
        
        function obj = setelev(obj, c, e)
            % SETELEV Assign elevation at coordinate.
            %   SETELEV(OBJ, C, E) assigns elevation E to the map tile at
            %   coordinate C.
            %
            %   C is an Mx2 matrix. Each row specifies an x-y coordinate.
            %   
            %   E is an M-element column vector. The map tile at coordinate
            %   C(m,:) will be assigned the elevation E(m).
            
            %% Validate input and output.
            % Check the input and output arguments.
            nargoutchk(0, 1)
            narginchk(3, 3)
            validateattributes(c, {'numeric'}, {'real','ncols',2}, '', 'C')
            validateattributes(e, {'numeric'}, ...
                {'real','ncols',1,'nrows',size(c,1)}, '', 'E')
            
            %% Set elevation at coordinate.
            i = obj.idx(c);
            if any(isnan(i))
                inan = i(find(isnan(i), 1));
                error(['Coordinate [', num2str(c(inan,1)), '; ', ...
                    num2str(c(inan,2)), '] lies outside map.'])
            end
            obj.elevation(i) = e;
        end
        
        function l = limits(obj)
            % LIMITS Spatial map extent.
            %   L = LIMITS(OBJ) returns a 2x2 matrix that specifies the
            %   extent of the elevation map in x- and y-direction:
            %      L = [xmin, xmax; ymin, ymax].
            
            nargoutchk(0, 1)
            narginchk(1, 1)
            s = repmat(obj.support, 2, 1);
            l = (s + [0,0; obj.extension*obj.resolution]).';
        end
        
        function m = mask(obj, p)
            % MASK Tell which tiles lie below or above given points.
            %   M = MASK(OBJ, P) returns a boolean IxJ matrix that 
            %   indicates whether an elevation map tile lies below or above
            %   one of the given points P or not. 
            %   I and J are the numbers of tiles in x- and y-direction, 
            %   respectively.
            %
            %   P is an Mx2 matrix, whose columns contain the x-y
            %   coordinates of the points, respectively.
            
            %% Validate input and output.
            nargoutchk(0, 1)
            narginchk(2, 2)
            validateattributes(p, {'numeric'}, {'real','ncols',2}, '', 'P')
            
            %% Compute mask.
            % Initialize the mask.
            m = false(size(obj.elevation));
            
            % Compute the indices of the tiles above or below which points
            % are located.
            i = idx(obj,p);
            
            % Set the mask values corresponding to these tiles to true.
            m(i(isfinite(i))) = true;
        end
        
        function d = diff(obj, p)
            % DIFF Height difference between 3D points and elevation map.
            %   D = DIFF(OBJ, P) computes the height difference between 
            %   elevation map OBJ and 3D points P.
            %
            %   P is an Mx3 matrix whose rows contain x-y-z coordinates.
            %
            %   D is an M-element column vector. D(m) yields the distance
            %   between P(m,:) and the corresponding map tile in 
            %   z direction. Positive values mean the point lies above the
            %   tile, negative values mean the point lies below the tile.
            %
            %   If there is no corresponding map tile for P(m,:), D(m)
            %   returns NaN.
            
            %% Validate input and output.
            nargoutchk(0, 1)
            narginchk(2, 2)
            validateattributes(p, {'numeric'}, {'real','ncols',3}, '', 'P')
            
            %% Compute difference in z.
            d = p(:,3) - obj.getelev(p(:,1:2));
        end
        
        function m = match(obj, pc)
            % MATCH Compute distance between point cloud and elevation map.
            %   M = MATCH(OBJ, PC) computes the mean distance M in 
            %   z-direction between the point cloud PC and the 
            %   elevation map OBJ.
            %
            %   PC is a point cloud object.
            %   M is a scalar that indicates the mean error in z-direction.
            %
            %   MATCH computes the z-distance between each point of PC and 
            %   the corresponding tile of the elevation map. 
            %   The computed distance is asymmetric: If the point lies 
            %   above the tile, the distance adds to the total error. If 
            %   the point lies below the tile, the distance is set to zero.
            %   This accounts for the fact that the map tiles represent the
            %   heighest points of the point cloud map only.
            
            %% Validate input.
            nargoutchk(0, 1)
            narginchk(2, 2)
            validateattributes(pc, {'pointCloud'}, {'scalar'}, '', 'PC')

            %% Compute mean distance.
            % Get the point coordinates as an Mx3 matrix.
            point = reshape(pc.Location, pc.Count, 3, 1);

            % Compute the mean of all positive distances between the points
            % and the elevation map.
            m = mean(constrain(diff(obj, point), [0, +Inf]));
        end
        
        function m = fillnan(obj, ws)
            % FILLNAN Estimate elevation values of NaN tiles.
            %   M = FILLNAN(OBJ) estimates the elevation values of NaN map
            %   tiles using a median filter.
            %
            %   M is an elevationmap object.
            %
            %   For each map tile of M, FILLNAN takes a window of 3x3 tiles
            %   centered around the corresponding tile of OBJ. 
            %   FILLNAN then computes the median of all tiles in the window
            %   and assigns the resulting value to the tile of M.
            %
            %   M = FILLNAN(OBJ, WS) performs median filtering using a
            %   user-defined window size WS.
            %
            %   WS is a 2-element vector containing odd integers. 
            %   The elements specify the extent of the window in x and y 
            %   direction, respectively.
            %
            %   Example:
            %      em = elevationmap(pcread('teapot.ply'), 0.05);
            %      em = em.fillnan([9,9])

            %% Validate input and output.
            % Check the number of input and output arguments.
            nargoutchk(0, 1)
            narginchk(1, 2)
            
            % If the window size is not defined, use the default.
            if nargin < 2
                ws = [3,3];
            end
            
            % Validate the window size.
            validateattributes(ws, {'numeric'}, ...
                {'odd', 'positive', 'numel',2}, '', 'WS')
            
            %% Filter values.
            % Compute half the window size.
            dn = floor(ws / 2);
            
            % Compute the size of the data matrix.
            sx = size(obj.elevation, 1);
            sy = size(obj.elevation, 2);

            % Compute an elevation matrix with its NaN values filled.
            % Loop over all map tiles.
            e = obj.elevation;
            parfor x = 1 : sx
                for y = 1 : sy
                    % If the map tile elevation is NaN, set it to the
                    % median of all tiles in the window.
                    if ~isfinite(e(x,y)) %#ok<*PFBNS>
                        % Compute the indices of the tiles in the window.
                        i = [max([1, x-dn(1); 1, y-dn(2)], [], 2), ...
                            min([sx, x+dn(1); sy, y+dn(2)], [], 2)]; 
                        
                        % Compute the median of the tiles in the window.
                        ew = obj.elevation(i(1,1):i(1,2),i(2,1):i(2,2));
                        e(x,y) = median(ew(:), 'omitnan');
                    end
                end
            end
            
            % Create the result elevation map.
            m = obj;
            m.elevation = e;
        end
        
        function m = setprior(obj, v)
            % SETPRIOR Set NaN tiles to value.
            %   M = SETPRIOR(OBJ, V) sets all NaN tiles to value V.
            
            %% Validate input and output.
            nargoutchk(0, 1)
            narginchk(2, 2)
            validateattributes(v, {'numeric'}, {'scalar', 'real'}, '', 'V')
            
            %% Set prior.
            m = obj;
            m.elevation(isnan(obj.elevation)) = v;
        end
        
        function plot(obj)
            % PLOT Visualize elevation map.
            
            %% Validate output.
            nargoutchk(0, 0)
            
            %% Plot.
            % Compute the x and y coordinates of each grid node.
            xgv = obj.support(1) ...
                + kron(0:obj.extension(1), [1,1]) * obj.resolution;
            ygv = obj.support(2) ...
                + kron(0:obj.extension(2), [1,1]) * obj.resolution;
            [x, y] = ndgrid(xgv(2:end-1), ygv(2:end-1));
            
            % Compute the elevation of each grid node.
            z = kron(obj.elevation, ones(2));
            
            % Plot the map surface.
            surf(x, y, z, 'EdgeColor', 'none')
            shading interp
            axis equal
            
            % Plot decoration.
            grid
            labelaxes
        end
    end 
end
