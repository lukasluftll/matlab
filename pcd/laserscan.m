classdef laserscan < handle
    % LASERSCAN Object for storing a 3D laser scan.
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS) creates a laser scan
    %   object.
    %
    %   AZIMUTH, ELEVATION, and RADIUS are HEIGHTxWIDTH matrices which
    %   define the rays the laser sensor measured in spherical coordinates.
    %   AZIMUTH and ELEVATION must be finite. RADIUS may contain infinite 
    %   values, which are interpreted as no-return rays.
    %
    %   The laser sensor pose is assumed to be equal to the reference 
    %   frame.
    %
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS, TFORM) additionally
    %   specifies the laser sensor pose with respect to the global 
    %   reference frame.
    %
    %   TFORM is a homogeneous 4x4 transformation matrix defining rotation 
    %   and translation only, no sheering, scaling, etc.
    %
    %   L = LASERSCAN(AZIMUTH, ELEVATION, RADIUS, TFORM, RLIM) additionally
    %   defines the minimum and maximum radius the sensor is able to 
    %   measure using the ordered 2-element vector RLIM.
    %
    %   Example:
    %      pcd = pcdread('castle.pcd');
    %      l = laserscan(pcd.azimuth, pcd.elevation, pcd.radius)
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
        
        function set.rlim(obj, rlim)
            if numel(rlim) ~= 2
                error('')
            end
            
            if rlim(1) < 0
                error('')
            end
            
            if diff(rlim) <= 0
                error('RLIM(2) must be greater than RLIM(1).');
            end
            
            obj.rlim = rlim;
        end
    end
    
    methods ( Access = public )
        function obj = laserscan(azimuth, elevation, radius, tform, rlim)
            narginchk(3, 5);
                      
            if nargin < 4
                tform = eye(4);
            end
            
            if nargin < 5
                rlim = [0; +Inf];
            end
            
            obj.azimuth = azimuth;
            obj.elevation = elevation;
            obj.radius = radius;
            obj.tform = tform;
            obj.rlim = rlim;
        end
        
        function inan = nan(obj)
            inan = ~isfinite(obj.radius) | obj.radius < obj.rlim(1) ...
                | obj.radius > obj.rlim(2);
        end
        
        function plot(obj)
            if isfinite(obj.rlim(2))
                rmax = max(obj.radius);
            else
                rmax = obj.rlim(2);
            end
            
            plotradius = obj.radius;
            plotradius(inan(obj)) = rmax;
            rayplot(obj.azimuth, obj.elevation, plotradius, inan(obj));
        end
    end
    
end

