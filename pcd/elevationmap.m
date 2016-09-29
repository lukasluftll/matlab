classdef elevationmap
    % ELEVATIONMAP Object for storing an elevation data grid.
    %   EM = ELEVATIONMAP(PC, RES) rasterizes the x-y plane with grid 
    %   resolution RES and sets the elevation of each grid tile to the 
    %   maximum z coordinate of all points of point cloud PC that belong
    %   to this tile.
    %
    %   PC is a pointCloud object. RES is a positive scalar.
    %
    %   Example:
    %      em = elevationmap(pcread('teapot.ply'), 0.1)
    %
    %   ELEVATIONMAP methods:
    %   GETELEV   - Elevation at coordinate
    %   LIMITS    - Extension of map
    %   MINUS     - Height difference between elevation map and 3D points
    %   MEDFILT   - 2D median filtering
    %   PLOT      - Visualize elevation map
    %
    %   See also POINTCLOUD.
    
    % Copyright 2016 Alexander Schaefer
    
    properties ( SetAccess = private )
        % SUPPORT Minimum x and y coordinates of the map; 1x2 vector.
        support
        
        % EXTENSION Extension of the map in numbers of map tiles in x and 
        % y direction; 1x2 vector.
        extension
        
        % RESOLUTION Edge length of each tile; scalar.
        resolution

        % ELEVATION Elevation data of the map; IxJ matrix with I and J 
        % being the number of tiles in x and y direction.
        elevation
    end
    
    methods ( Access = private )
        function i = idx(obj, p)
            % IDX Compute tile index given x-y coordinates.
            %   I = IDX(OBJ, P) computes the linear index I of the tile to
            %   which the given x-y coordinates P belong.
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
            d = p(valid,:) - repmat(obj.support, sum(valid), 1);
            isub = floor(d / obj.resolution) + 1;
            i(valid) = sub2ind(obj.extension, isub(:,1), isub(:,2));
            
        end                
    end
    
    methods ( Access = public )      
        function obj = elevationmap(pc, res)
            % ELEVATIONMAP Constructor.
            
            %% Validate input and output.
            % Check number of input and output arguments.
            nargoutchk(0, 1)
            narginchk(2, 2)
            
            % Check the validity of the point cloud.
            validateattributes(pc, {'pointCloud'}, {'numel',1}, '', 'PC')            
            validateattributes(pc.XLimits, {'numeric'}, {'finite'}, ...
                '', 'PC.XLIMITS')
            validateattributes(pc.YLimits, {'numeric'}, {'finite'}, ...
                '', 'PC.YLIMITS')
            
            % Check the validity of the resolution.
            validateattributes(res, {'numeric'}, ...
                {'scalar', 'nonnegative'}, '', 'RES')

            %% Compute property values.
            % Store the resolution.
            obj.resolution = res;
            
            % Compute the minimum x and y coordinates of the map.
            obj.support = floor([pc.XLimits(1),pc.YLimits(1)] / res) * res;
            
            % Compute the number of map tiles in x and y direction.
            obj.extension = floor(...
                ([pc.XLimits(2),pc.YLimits(2)]-obj.support) / res) + 1;
            
            % Set the elevation data to the minimum elevation.
            obj.elevation = NaN(obj.extension);
            
            % Remove all NaN and infinite points from the point cloud.
            pc = removeInvalidPoints(pc);
            
            % Make sure the point cloud is represented by an Mx3 vector.
            point = reshape(pc.Location(:), pc.Count, 3, 1);
            
            % Compute the indices of the tiles to which each point belongs.
            ie = obj.idx(point(:,1:2));
            
            % Set the elevation of each tile to the maximum z coordinate of
            % all points belonging to this tile.
            for i = 1 : size(point,1)
                obj.elevation(ie(i))=max(obj.elevation(ie(i)),point(i,3));
            end
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
            
            % Look up the elevation for each coordinate. 
            i = obj.idx(xy);
            fin = isfinite(i);
            e(fin) = obj.elevation(i(fin));
        end
        
        function l = limits(obj)
            % LIMITS Map extent.
            %   L = LIMITS(OBJ) returns a 2x2 matrix that specifies the
            %   extent of the elevation map:
            %      L = [xmin, xmax; ymin, ymax].
            
            nargoutchk(0, 1)
            s = repmat(obj.support, 2, 1);
            l = (s + [0,0; obj.extension*obj.resolution]).';
        end
        
        function d = minus(obj, p)
            % -  Height difference between elevation map and 3D points.
            %   D = OBJ - P computes the height difference between 
            %   elevation map OBJ and 3D points P.
            %
            %   P is an Mx3 matrix whose rows contain x-y-z coordinates.
            %
            %   D is an M-element column vector. D(m) yields the height
            %   difference between P(m,:) and the corresponding map tile.
            %   If there is no corresponding map tile for P(m,:), D(m)
            %   returns NaN.
            
            %% Validate input and output.
            nargoutchk(0, 1)
            narginchk(2, 2)
            validateattributes(p, {'numeric'}, {'real','ncols',3}, '', 'P')
            
            %% Compute height difference.
            d = obj.z(p(:,1:2)) - p(:,3);
        end
        
        function m = medfilt(obj, ws)
            % MEDFILT 2D median filtering.
            %   M = MEDFILT(OBJ) performs median filtering of the
            %   elevation map OBJ using a window of 3x3 elements. The
            %   window is always centered around the cell whose value is
            %   computed.
            %
            %   M is an elevationmap object. The elevation value of each
            %   cell is equal to the median of the values of OBJ inside the
            %   window.
            %
            %   M = MEDFILT(OBJ, WS) performs median filtering using a
            %   user-defined window size WS.
            %
            %   WS is a 2-element vector containing odd integers. 
            %   The elements specify the extent of the window in x and y 
            %   direction, respectively.

            %% Validate input and output.
            % Check the numer of input and output arguments.
            nargoutchk(0, 1)
            narginchk(1, 2)
            
            % If the window size is not defined, do it now.
            if nargin < 2
                ws = [3,3];
            end
            
            % Validate the window size.
            validateattributes(ws, {'numeric'}, ...
                {'odd', 'positive', 'numel',2}, '', 'WS')
            
            %% Perform filtering.
            % Compute the number of neighbors in each direction.
            nn = (ws-1) / 2;
            m = obj;
            m.elevation = median(obj.elevation(
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
