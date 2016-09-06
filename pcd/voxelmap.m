classdef voxelmap
    % VOXELMAP Object for storing a 3D voxel map.
    %   VM = VOXELMAP(DATA) creates a voxel map object.
    %   
    %   DATA is an IxJxK matrix that contains the map data.
    %
    %   VM = VOXELMAP(DATA, XGV, YGV, ZGV) creates a voxel map object with
    %   a defined voxel grid.
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
    %   VM = VOXELMAP(DATA, XGV, YGV, ZGV, PRIOR) defines the value of all
    %   points outside the volume covered by the voxelmap. The default 
    %   value for PRIOR is NaN.
    %
    %   Example:
    %      vm = voxelmap(rand(10,10,10))
    %
    %   VOXELMAP properties:
    %   DATA     - Map data
    %   XGV      - x-axis grid vector
    %   YGV      - y-axis grid vector
    %   ZGV      - z-axis grid vector
    %   PRIOR    - Value of all points outside the map volume
    %
    %   VOXELMAP methods:
    %   PLUS      - Add map data element-wise
    %   SUM       - Sum up map data element-wise
    %   MINUS     - Subtract map data element-wise
    %   TIMES     - Multiply map data element-wise
    %   RDIVIDE   - Divide map data element-wise from right
    %   LOG       - Natural logarithm of map dta
    %   UMINUS    - Unary minus of map data
    %   CONSTRAIN - Fit map data to interval
    %   SUBVOLUME - Map subvolume
    %   LIMITS    - Spatial limits of map 
    %   HISTOGRAM - Histogram of map data
    %   PLOT      - Visualize map using transparency values
    %
    %   See also LFMAP, REFMAP, DECAYMAP.
    
    % Copyright 2016 Alexander Schaefer
    
    properties 
        data;
        prior;
    end
    
    properties ( SetAccess = private )
        xgv;
        ygv;
        zgv;
    end
    
    methods
        function obj = set.data(obj, data)
            % SET.DATA Change map data.
            
            % Check if DATA has the same size as the current map data.
            if ~isempty(obj.data) && (ndims(data) ~= ndims(obj.data) ...
                    || any(size(data) ~= size(obj.data)))
                error(['DATA must be a matrix of size ', ...
                    num2str(size(obj.data, 1)), 'x', ...
                    num2str(size(obj.data, 2)), 'x', ...
                    num2str(size(obj.data, 3)), '.'])
            end
            
            % Assign new data to object.
            obj.data = data;
        end
        
        function obj = set.prior(obj, prior)
            % SET.PRIOR Set map prior.
            
            % Check the size and type of PRIOR.
            if numel(prior) ~= 1 || ~isnumeric(prior)
                error('PRIOR must be a numeric scalar.')
            end
            
            % Assign new prior to object.
            obj.prior = prior;
        end
        
        function chkewo(a, b)
            % CHKEWO Check if two voxelmaps are amenable to element-wise 
            % operations.
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Check whether both input arguments are voxelmap objects.
            if ~isa(a, 'voxelmap') || ~isa(b, 'voxelmap')
                error('A and B must be voxelmap objects.')
            end
            
            %% Check readiness for element-wise operations.
            % Check if the grid vectors match.
            if ~(all(a.xgv == b.xgv) && all(a.ygv == b.ygv) ...
                && any(a.zgv == b.zgv))
                error('Grid vectors of both voxelmaps must match.')
            end
                        
            % Check if the data matrix sizes match.
            if ndims(a.data) ~= ndims(b.data) ...
                    || ~all(size(a.data) == size(b.data))
                error('Voxelmaps must be of same size.');
            end
        end
        
        function [a, b] = any2voxelmap(a, b)
            % ANY2VOXELMAP Convert to voxelmap.
            %   [A, B] = ANY2VOXELMAP(A, B) returns two voxelmap objects A 
            %   and B. At least one of the input arguments must already be
            %   a voxelmap object.

            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Check whether at least on input argument is a voxelmap
            % object.
            if ~isa(a, 'voxelmap') && ~isa(b, 'voxelmap')
                error('At least one input must be a voxelmap object.')
            end
            
            %% Convert.
            if ~isa(a, 'voxelmap')
                a = voxelmap(a, b.xgv, b.ygv, b.zgv, b.prior);
            elseif ~isa(b, 'voxelmap')
                b = voxelmap(b, a.xgv, a.ygv, a.zgv, a.prior);
            end
        end
        
        function [a, b] = fitsize(a, b)
            % FITSIZE Fit size of map data to match, if possible.
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Make sure both inputs are voxelmap objects.
            if ~isa(a, 'voxelmap') || ~isa(b, 'voxelmap')
                error('A and B must be voxelmap objects.')
            end
            
            %% Fit size.
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            if isempty(a.data)
                a.data = NaN(size(b.data));
            elseif isempty(b.data)
                b.data = NaN(size(a.data));
            end
        end
    end
        
    methods ( Access = public )
        function obj = voxelmap(data, xgv, ygv, zgv, prior)
            % VOXELMAP Constructor.
            
            %% Validate input.
            % Check number of input arguments.
            switch nargin
                case 1
                    xgv = 0 : size(data, 1);
                    ygv = 0 : size(data, 2);
                    zgv = 0 : size(data, 3);
                    prior = NaN;
                case 4
                    prior = NaN;
                case 5
                otherwise
                    error('VOXELMAP takes 1, 4, or 5 input arguments.')
            end
            
            % Check the grid vectors.
            gvchk(xgv, ygv, zgv);
            
            %% Assign input to properties.
            obj.data = data;
            obj.prior = prior;
            obj.xgv = xgv;
            obj.ygv = ygv;
            obj.zgv = zgv;
        end
        
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
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Convert the input arguments to voxelmap.
            [a, b] = any2voxelmap(a, b);
            
            % If any voxelmap is empty, set its size to match the other
            % voxelmap.
            [a, b] = fitsize(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            %% Add map data.
            c = voxelmap(a.data + b.data, a.xgv, a.ygv, a.zgv, ...
                a.prior + b.prior);
        end
        
        function b = sum(a)
            % SUM Sum of map data of multiple voxelmaps.
            %   S = SUM(A) returns a voxelmap object whose map data
            %   contains the element-wise sum of the map data matrices in 
            %   the voxelmap array A.
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(1, 1)
            
            % Check whether element-wise operations can be performed.
            for i = 2 : numel(a)
                chkewo(a(1), a(i))
            end
            
            %% Compute sum.
            b = a(1);
            for i = 2 : numel(a)
                b = b + a(i);
            end
        end
            
        function c = minus(a, b)
            % -  Subtract map data element-wise.
            %   C = A - B subtracts the map data B element-wise from A and 
            %   returns the difference.
            %
            %   A and B are voxelmap objects. One of them may also be a 3D
            %   matrix of the same size as the voxelmap object's data
            %   matrix.
            %
            %   C is a voxelmap object.
            c = plus(a, -b);
        end
        
        function c = times(a, b)
            % .*  Multiply map data element-wise.
            %   C = A .* B multiplies the map data A element-wise by B and 
            %   returns the product C.
            %
            %   A and B are voxelmap objects. One of them may also be a 3D
            %   matrix of the same size as the voxelmap object's data
            %   matrix.
            %
            %   C is a voxelmap object. 
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
                        
            % Force the input argument to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            %% Multiply map data.
            c = voxelmap(a.data .* b.data, a.xgv, a.ygv, a.zgv, ...
                a.prior * b.prior);
        end
        
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
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
                        
            % Force the input argument to be voxelmaps.
            [a, b] = any2voxelmap(a, b);
            
            % Check whether element-wise operations can be performed.
            chkewo(a, b)
            
            %% Divide cell data.
            c = voxelmap(a.data ./ b.data, a.xgv, a.ygv, a.zgv, ...
                a.prior / b.prior);
        end
        
        function obj = log(obj)
            % LOG(OBJ) Natural logarithm.
            %   L = LOG(OBJ) computes the natural logarithm for every voxel.
            
            % Check number of input arguments.
            narginchk(1, 1)
            
            % Compute logarithm.
            obj.data = log(obj.data);
            obj.prior = log(obj.prior);
        end
        
        function obj = uminus(obj)
            % -  Unary minus.
            %   Negates the values of all map elements of OBJ.
            
            % Check number of input arguments.
            narginchk(1, 1)
            
            % Negate data.
            obj.data = -obj.data;
            obj.prior = -obj.prior;
        end
        
        function obj = constrain(obj, lim)
            % CONSTRAIN Fit map data to interval.
            %   OBJ = CONSTRAIN(OBJ, LIM) fits all map data values of
            %   voxelmap object OBJ into the interval defined by the
            %   ordered 2-element vector LIM.
            %
            %   NaN values are set to LIM(1).
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
            
            %% Fit data to interval.
            obj.data = constrain(obj.data, lim);
            obj.data(isnan(obj.data)) = lim(1);
            obj.prior = constrain(obj.prior, lim);
            obj.prior(isnan(obj.prior)) = lim(1);
        end
        
        function submap = subvolume(obj, lim)
            % SUBVOLUME Get part of given map.
            %   OBJ = SUBVOLUME(OBJ, LIM) extracts a subset of the 
            %   voxelmap object OBJ using the specified axis-aligned limits
            %   LIM. 
            %
            %   LIM is a vector [xmin, xmax, ymin, ymax, zmin, zmax].
            %   Any NaN values in the limits indicate that the voxelmap 
            %   should not be cropped along that axis.
            %   The returned subvolume of the map consists of all voxels
            %   that lie completely within LIM.
            
            %% Validate input.
            % Check number of input arguments.
            narginchk(2, 2)
            
            % Check the limits.
            validateattributes(lim,{'numeric'},{'numel',6},'','LIM')
            
            % Set all NaN values to minimum or maximum limits.
            mapLimits = obj.limits;
            lim(isnan(lim)) = mapLimits(isnan(lim));
            
            %% Create map subvolume.
            % Compute the logical indices of the grid vectors of the
            % subvolume.
            ix = lim(1)<=obj.xgv & obj.xgv<=lim(2);
            iy = lim(3)<=obj.ygv & obj.ygv<=lim(4);
            iz = lim(5)<=obj.zgv & obj.zgv<=lim(6);
            
            % Compute grid vectors of the subvolume.
            subxgv = obj.xgv(ix);
            subygv = obj.ygv(iy);
            subzgv = obj.zgv(iz);
            
            % Compute the logical indices of the voxels in the subvolume.
            ix(find(ix, 1, 'last')) = false;
            iy(find(iy, 1, 'last')) = false;
            iz(find(iz, 1, 'last')) = false;
            
            % Get the subvolume of the map data.
            [ix,iy,iz] = ndgrid(ix,iy,iz);
            subgridsize = [numel(subxgv),numel(subygv),numel(subzgv)] - 1;
            subdata = reshape(obj.data(ix&iy&iz), subgridsize);
            
            % Create the subvolume voxelmap.
            submap = voxelmap(subdata, subxgv, subygv, subzgv, obj.prior);
        end
        
        function l = limits(obj)
            % LIMITS Spatial limits of map.
            %   L = LIMITS(OBJ) returns a 6-element row vector that
            %   contains the spatial limits of the voxelmap object OBJ in
            %   the order [xmin, xmax, ymin, ymax, zmin, zmax].
            
            l = [min(obj.xgv), max(obj.xgv), ...
                min(obj.ygv), max(obj.ygv), min(obj.zgv),max(obj.zgv)];
        end
        
        function h = histogram(obj, varargin)
            % HISTOGRAM Plot histogram of map data.
            %   HISTOGRAM(OBJ) plots a histogram of all map data contained
            %   in the voxelmap object OBJ. The function works like the 
            %   standard HISTORGRAM function.
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

        function plot(obj, mode)
            % PLOT Plot voxel map.
            %   PLOT(OBJ) visualizes the voxelmap object OBJ using
            %   semi-transparent voxels. 
            %
            %   PLOT(OBJ, MODE) determines how the map data is represented
            %   as colors. MODE can assume the following values:
            %
            %      'scale' (default) - The highest value of OBJ is drawn as
            %                          completely intransparent voxel, the 
            %                          lowest value of OBJ is drawn as
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
            
            %% Convert data to colors.
            % Compute minimum and maximum data elements.
            lim = [min(plotdata(:)), max(plotdata(:))];
            
            % Compute transparency values to plot
            switch lower(mode)
                case 'scale'
                    plotdata = (plotdata - lim(1)) / diff(lim);
                case 'direct'
                    if lim(1) < 0 || lim(2) > 1
                        warning('Map data exceeds [0;1].')
                    end
                otherwise
                    error(['MODE ', mode, ' not supported.'])
            end
            
            % Crop data to interval [0; 1].
            plotdata = constrain(plotdata, [0,1]);
            
            %% Plot.
            alphaplot(plotdata, plotdata, obj.xgv, obj.ygv, obj.zgv);
            axis image vis3d
            grid on
            labelaxes
        end
    end
end
