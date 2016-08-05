classdef laserscan < handle
    
    properties
        tform;
        azimuth;
        elevation;
        radius;
        rlim;
    end
    
    methods
        function set.azimuth(obj, azimuth)
            if ~isempty(obj.azimuth) && size(azimuth) ~= size(obj.azimuth)
                error('')
            end
            
            if ~all(isfinite(azimuth(:)))
                error('')
            end
            
            obj.azimuth = azimuth;
        end
        
        function set.elevation(obj, elevation)
            if ~isempty(obj.elevation) ...
                    && size(elevation)~=size(obj.elevation)
                error('')
            end
            
            if ~all(isfinite(elevation(:)))
                error('')
            end
            
            obj.elevation = elevation;
        end
        
        function set.radius(obj, radius)
            if ~isempty(obj.radius) && size(radius) ~= size(obj.radius)
                error('')
            end
            
            obj.radius = radius;
        end
        
        function set.tform(obj, tform)
            if ~ishtform(tform)
                error('')
            end
            
            obj.tform = tform;
        end
    end
    
    methods ( Access = public )
        function obj = laserscan(azimuth, elevation, radius, tform)
            narginchk(3, 4);
                      
            if nargin < 4
                tform = eye(4);
            end
            
            obj.azimuth = azimuth;
            obj.elevation = elevation;
            obj.radius = radius;
            obj.tform = tform;
        end
    end
    
end

