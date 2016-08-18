classdef voxelmap < handle
    % VOXELMAP Object for storing a 3D voxel map.
    %   M = VOXELMAP(DATA) creates a voxel map object.
    %   
    %   DATA is a IxJxK matrix that contains the map data.
    %
    %   M = VOXELMAP(DATA, XGV, YGV, ZGV) creates a voxel map object with a
    %   defined voxel grid.
    %   
    %   XGV, YGV, ZGV are vectors that define the rasterization of the 
    %   voxel grid. A voxel with index [i,j,k] contains all points [x,y,z]
    %   that satisfy the inequality:
    %
    %      (XGV(i) <= x < XGV(i+1))
    %      && (YGV(j) <= y < YGV(j+1)) 
    %      && (ZGV(k) <= z < ZGV(k+1))
    %
    %   DATA is a IxJxK matrix that contains the map data, where 
    %   I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
    %
    %   Example:
    %      data = [0.1, 0.2, 0.3; 0, 0, 0; 0, 0, 0];
    %      data = cat(3, data, zeros(3), 0.5 * ones(3));
    %      m = voxelmap(data)
    %
    %   VOXELMAP properties:
    %   DATA     - Cell data
    %   XGV      - x-axis grid vector
    %   YGV      - y-axis grid vector
    %   ZGV      - z-axis grid vector
    %
    %   VOXELMAP methods:
    %   PLUS     - Add cell data element-wise
    %   MINUS    - Subtract cell data element-wise
    %   TIMES    - Multiply cell data element-wise
    %   RDIVIDE  - Divide cell data element-wise from right
    %   PLOT     - Visualize the map using transparency values
    %
    %   See also LFMAP, REFMAP, DECAYMAP.
    
    % Copyright 2016 Alexander Schaefer
    
    properties 
        data;
    end
    
    properties ( SetAccess = private )
        xgv;
        ygv;
        zgv;
    end
    
    methods
        % Change the map data.
        function set.data(obj, data)
            % Check if DATA has the same size as the current map data.
            if ~isempty(obj.data) && any(size(data) ~= size(obj.data))
                error(['DATA must be a matrix of size ', ...
                    num2str(size(obj.data, 1)), 'x', ...
                    num2str(size(obj.data, 2)), 'x', ...
                    num2str(size(obj.data, 3)), '.'])
            end
            
            % Assign new map data to object.
            obj.data = data;
        end
        
        % Check if the data matrices of the two voxelmaps are amenable to 
        % element-wise operations.
        function chkewo(l, r)
            % Check if the grid vectors match.
            if ~(all(l.xgv == r.xgv) && all(l.ygv == r.ygv) ...
                && any(l.zgv == r.zgv))
                error('Grid vectors of both voxelmaps must match.')
            end
                        
            % Check if the data matrix sizes match.
            if ~all(size(l.data) == size(r.data))
                error('Voxelmaps must be of same size.');
            end
        end
        
        % Converts a pair of input arguments to voxelmap objects.
        % At least one of the arguments must already be a voxelmap object.
        function [a, b] = any2voxelmap(a, b)
            % If any input argument is not a voxelmap, convert it.
            if ~isa(a, 'voxelmap')
                a = voxelmap(a, b.xgv, b.ygv, b.zgv);
            elseif ~isa(b, 'voxelmap')
                b = voxelmap(b, a.xgv, a.ygv, a.zgv);
            end
        end
    end
        
    methods ( Access = public )
        % Construct a voxelmap object.
        function obj = voxelmap(data, xgv, ygv, zgv)
            % Check number of input arguments.
            if ~(nargin == 1 || nargin == 4)
                error('VOXELMAP takes 1 or 4 input arguments.')
            end
            
            % If the grid vectors are not given, construct them.
            if nargin == 1
                xgv = 0 : size(data, 1);
                ygv = 0 : size(data, 2);
                zgv = 0 : size(data, 3);
            end
            
            % Check the grid vectors.
            gvchk(xgv, ygv, zgv);
            
            % Store the input.
            obj.data = data;
            obj.xgv = xgv;
            obj.ygv = ygv;
            obj.zgv = zgv;
        end
        
        % Add cell data element-wise.
        function c = plus(a, b)
            % +  Add cell data element-wise.
            
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Force the input arguments to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            if isempty(a.data)
                a.data = zeros(size(b.data));
            elseif isempty(b.data)
                b.data = zeros(size(a.data));
            end
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            % Add cell data.
            c = voxelmap(a.data + b.data, a.xgv, a.ygv, a.zgv);
        end
            
        % Subtract cell data element-wise.
        function c = minus(a, b)
            % -  Subtract cell data element-wise.
            
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Force the input argument to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            if isempty(a.data)
                a.data = zeros(size(b.data));
            elseif isempty(b.data)
                b.data = zeros(size(a.data));
            end
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            % Subtract cell data.
            c = voxelmap(a.data - b.data, a.xgv, a.ygv, a.zgv);
        end
        
        % Multiply cell data element-wise.
        function c = times(a, b)
            % .*  Multiply cell data element-wise.
            
            % Check number of input arguments.
            narginchk(2, 2)
                        
            % Force the input argument to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            % Multiply cell data.
            c = voxelmap(a.data .* b.data, a.xgv, a.ygv, a.zgv);
        end
        
        % Divide cell data element-wise from right.
        function c = rdivide(a, b)
            % ./  Divide cell data element-wise from right.
            
            % Check number of input arguments.
            narginchk(2, 2)
                        
            % Force the input argument to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            % Divide cell data.
            c = voxelmap(a.data ./ b.data, a.xgv, a.ygv, a.zgv);
        end

        % Plot the map.
        function plot(obj, mode)
            % PLOT Plot voxel map.
            %   PLOT(OBJ, MODE) visualizes the voxelmap object VM using
            %   semi-transparent voxels. 
            %
            %   MODE can assume the following values:
            %
            %      'scale' (default) - The highest value of VM is drawn as
            %                          completely intransparent voxel, the 
            %                          lowest value of VM is drawn as
            %                          completely transparent voxel.
            %      'direct'          - Data values are directly interpreted 
            %                          as transparencies, with 0
            %                          representing complete intransparency
            %                          and 1 representing complete
            %                          transparency.
            
            %% Validate input.
            % Check the number of input arguments.
            narginchk(1, 2);
            
            % If there is mode given, set it.
            if nargin < 2
                mode = 'scale';
            end
            
            % Set all infinite values to NaN.
            plotdata = obj.data;
            plotdata(~isfinite(plotdata)) = NaN;
            
            %% Plot map.
            switch lower(mode)
                case 'scale'
                    % Compute minimum and maximum data elements.
                    lim = [min(plotdata(:)), max(plotdata(:))];
                    
                    % Compute the transparency values to plot.
                    plotdata = (plotdata + lim(1)) / diff(lim);
                case 'direct'
                otherwise
                    error(['MODE ', mode, ' not supported.'])
            end
            
            % Crop data to interval [0; 1].
            plotdata = constrain(plotdata, [0,1]);
            
            % Plot.
            alphaplot(plotdata, obj.xgv, obj.ygv, obj.zgv);
            axis equal
            grid on
            xlabel('x')
            ylabel('y')
            zlabel('z')
        end
    end
end
