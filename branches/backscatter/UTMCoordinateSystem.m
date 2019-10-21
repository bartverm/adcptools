classdef UTMCoordinateSystem < ProjectedCoordinateSystem
    properties (Constant)
        description = 'UTM';
    end
    properties
        zone
    end
    methods
        function [x, y]=xy(obj,lat, lon)
            if isempty(obj.zone)
                [x,y,obj.zone]=geo2utm(lat,lon);
            else
                [x,y]=geo2utm(lat,lon,obj.zone);
            end
        end
    end
    
end