classdef gridmap
    
    properties
        data;
        support;
        edgelength;
        extent;
        prior;
    end
    
    methods (Access = protected)
        function i = idx(obj, p)
            % IDX Compute tile index given coordinates.
            %   I = IDX(OBJ, P) computes the linear index I of the tile to
            %   which the given x-y coordinates P correspond.
            %
            %   P is an MxD vector, where D is the number of dimensions of
            %   OBJ.DATA.
            %
            %   I is an M-element column vector. If P(m,:) lies outside the
            %   map, I(m) is NaN.
            
            %% Validate input and output.
            nargoutchk(0, 1)
            narginchk(2, 2)

            %% Compute tile index.
            d = p - repmat(obj.support, size(p,1), 1);
            isub = floor(d / obj.edgelength) + 1;
            isub(isub < 1 | isub > numel(obj.data)) = NaN;
            isub = num2cell(isub,1);
            i = sub2ind(obj.extent, isub{:});
        end                
    end
    
    methods
        function setdata(obj, d, p)
        end
        
        function d = getdata(obj, p)
        end
        
        function n = ndims(obj)
        end
        
        function l = limits(obj)
        end
        
        function y = plus(x1, x2)
        end
        
        function y = sum(x)
        end
        
        function y = minus(x1, x2)
        end
        
        function y = times(x1, x2)
        end
        
        function y = rdivide(x1, x2)
        end
        
        function obj = log(obj)
        end
        
        function obj = uminus(obj)
        end
        
        function obj = constrain(obj, lim)
        end
        
        function h = histogram(obj, varargin)
        end
    end
    
end

