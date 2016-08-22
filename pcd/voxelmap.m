classdef voxelmap < matlab.mixin.Copyable
    % VOXELMAP Object for storing a 3D voxel map.
    %   M = VOXELMAP(DATA) creates a voxel map object.
    %   
    %   DATA is an IxJxK matrix that contains the map data.
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
    %   DATA is an IxJxK matrix that contains the map data, where 
    %   I = numel(XGV)-1, J = numel(YGV)-1, and K = numel(ZGV)-1.
    %
    %   Example:
    %      data = rand(10, 10, 10);
    %      m = voxelmap(data)
    %
    %   VOXELMAP properties:
    %   DATA     - Map data
    %   XGV      - x-axis grid vector
    %   YGV      - y-axis grid vector
    %   ZGV      - z-axis grid vector
    %
    %   VOXELMAP methods:
    %   PLUS      - Add map data element-wise
    %   ADD       - Increment map data
    %   SUM       - Sum up map data element-wise
    %   MINUS     - Subtract map data element-wise
    %   SUBTRACT  - Decrement map data
    %   TIMES     - Multiply map data element-wise
    %   RDIVIDE   - Divide map data element-wise from right
    %   LOG       - Natural logarithm of map dta
    %   UMINUS    - Unary minus of map data
    %   CONSTRAIN - Fit map data to interval
    %   HISTOGRAM - Histogram of map data
    %   PLOT      - Visualize map using transparency values
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
        
        % Check if the map data of two voxelmaps are amenable to 
        % element-wise operations.
        function chkewo(a, b)
            % Check if the grid vectors match.
            if ~(all(a.xgv == b.xgv) && all(a.ygv == b.ygv) ...
                && any(a.zgv == b.zgv))
                error('Grid vectors of both voxelmaps must match.')
            end
                        
            % Check if the data matrix sizes match.
            if ~all(size(a.data) == size(b.data))
                error('Voxelmaps must be of same size.');
            end
        end
        
        % Convert a pair of input arguments to voxelmap objects.
        % At least one of the arguments must already be a voxelmap object.
        function [a, b] = any2voxelmap(a, b)
            % If any input argument is not a voxelmap, convert it.
            if ~isa(a, 'voxelmap')
                a = voxelmap(a, b.xgv, b.ygv, b.zgv);
            elseif ~isa(b, 'voxelmap')
                b = voxelmap(b, a.xgv, a.ygv, a.zgv);
            end
        end
        
        % Fit the size of the map data of two voxelmaps to match, if 
        % possible.
        function [a, b] = fitsize(a, b)
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            if isempty(a.data)
                a.data = zeros(size(b.data));
            elseif isempty(b.data)
                b.data = zeros(size(a.data));
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
        
        % Add map data element-wise.
        function c = plus(a, b)
            % +  Add map data element-wise.
            %   C = A + B adds the map data A and B element-wise and 
            %   returns the sum.
            %
            %   A and B are voxelmap objects. One of them may also be a 3D
            %   matrix of the same size as the voxelmap object's data
            %   matrix.
            %
            %   C is a voxelmap object.
            
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Convert the input arguments to voxelmap.
            [a, b] = any2voxelmap(a, b);
            
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            [a, b] = fitsize(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            % Add map data and create resulting voxelmap.
            c = voxelmap(a.data + b.data, a.xgv, a.ygv, a.zgv);
        end
        
        % Increment map data element-wise.
        function add(obj, d)
            % ADD Increment map data element-wise.
            %   ADD(OBJ, D) adds the data D to the map data of voxelmap
            %   object OBJ.
            %
            %   OBJ must be a voxelmap object.
            %
            %   D is a voxelmap object or a 3D matrix. The data matrix must
            %   be of the same size as the map data of OBJ.
            
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Convert the data to add to a voxelmap object.
            [obj, d] = any2voxelmap(obj, d);
            
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            [obj, d] = fitsize(obj, d);
            
            % Check whether element-wise operation can be performed.
            chkewo(obj, d)
            
            % Add the data to the voxelmap object.
            obj.data = obj.data + d.data;
        end
        
        % Sum of multiple voxelmaps.
        function s = sum(x)
            % SUM Sum of map data of multiple voxelmaps.
            %   S = SUM(X) returns a voxelmap object whose map data
            %   contains the element-wise sum of the map data matrices in 
            %   the voxelmap array X.
            
            % Check number of input arguments.
            narginchk(1, 1)
            
            % Check whether element-wise operations can be performed.
            for i = 2 : numel(x)
                chkewo(x(1), x(i))
            end
            
            % Compute the sum of the map data.
            s = x(1);
            for i = 2 : numel(x)
                s.add(x(i));
            end
        end
            
        % Subtract cell data element-wise.
        function c = minus(a, b)
            % -  Subtract cell data element-wise.
            %   C = A - B subtracts the map data B element-wise from A and 
            %   returns the difference.
            %
            %   A and B are voxelmap objects. One of them may also be a 3D
            %   matrix of the same size as the voxelmap object's data
            %   matrix.
            %
            %   C is a voxelmap object.   
            
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Force the input argument to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            [a, b] = fitsize(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            % Subtract cell data.
            c = voxelmap(a.data - b.data, a.xgv, a.ygv, a.zgv);
        end
        
        % Decrement map data element-wise.
        function subtract(obj, d)
            % SUBTRACT Decrement map data element-wise.
            %   SUBTRACT(OBJ, D) subtracts the data D from the map data of 
            %   voxelmap object OBJ.
            %
            %   OBJ must be a voxelmap object.
            %
            %   D is a voxelmap object or a 3D matrix. The data matrix must
            %   be of the same size as the map data of OBJ.
            add(obj, -d);
        end
        
        % Multiply cell data element-wise.
        function c = times(a, b)
            % .*  Multiply cell data element-wise.
            %   C = A .* B multiplies the map data A element-wise by B and 
            %   returns the product C.
            %
            %   A and B are voxelmap objects. One of them may also be a 3D
            %   matrix of the same size as the voxelmap object's data
            %   matrix.
            %
            %   C is a voxelmap object. 
            
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
            %   C = A ./ B divides the map data A element-wise from right 
            %   by B and returns the quotient.
            %
            %   A and B are voxelmap objects. One of them may also be a 3D
            %   matrix of the same size as the voxelmap object's data
            %   matrix.
            %
            %   C is a voxelmap object. 
            
            % Check number of input arguments.
            narginchk(2, 2)
                        
            % Force the input argument to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            % Divide cell data.
            c = voxelmap(a.data ./ b.data, a.xgv, a.ygv, a.zgv);
        end
        
        % Natural logarithm.
        function l = log(obj)
            % LOG(OBJ) Natural logarithm.
            %   L = LOG(OBJ) computes the natural logarithm of every cell.
            
            % Check number of input arguments.
            narginchk(1, 1)
            
            % Compute logarithm.
            l = voxelmap(log(obj.data), obj.xgv, obj.ygv, obj.zgv);
        end
        
        % Unary minus.
        function u = uminus(obj)
            % -  Unary minus.
            %   Negates the values of all map elements of OBJ.
            
            % Check number of input arguments.
            narginchk(1, 1)
            
            % Negate data.
            u = voxelmap(-obj.data, obj.xgv, obj.ygv, obj.zgv);
        end
        
        % Fit map data to interval.
        function y = constrain(obj, lim)
            % CONSTRAIN Fit map data to interval.
            %   Y = CONSTRAIN(OBJ, LIM) fits all map data values of
            %   voxelmap object OBJ into the interval defined by the
            %   ordered 2-element vector LIM.
            %
            %   NaN values are set to LIM(1).
            y = copy(obj);
            y.data = constrain(y.data, lim);
            y.data(isnan(y.data)) = lim(1);
        end
        
        % Plot histogram of map data.
        function h = histogram(obj, varargin)
            % HISTOGRAM Plot histogram of map data.
            %   HISTOGRAM(OBJ) plots a histogram of all map data contained
            %   in the voxelmap object OBJ. Apart from its first input
            %   argument, the function works exactly like the standard
            %   HISTORGRAM function.
            %
            %   Example:
            %      histogram(voxelmap(rand(10,10,10)))
            %
            % See also HISTOGRAM.
            if nargout > 0
                h = histogram(obj.data(:), varargin{:});
            else
                histogram(obj.data(:), varargin{:});
            end
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
                    plotdata = (plotdata - lim(1)) / diff(lim);
                case 'direct'
                otherwise
                    error(['MODE ', mode, ' not supported.'])
            end
            
            % Crop data to interval [0; 1].
            plotdata = constrain(plotdata, [0,1]);
            
            % Plot.
            alphaplot(plotdata, plotdata, obj.xgv, obj.ygv, obj.zgv);
            colormap('jet')
            axis image vis3d
            grid on
            labelaxes
        end
    end
end
