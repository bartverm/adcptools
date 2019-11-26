classdef UTMCoordinateSystem < ProjectedCoordinateSystem
% UTM coordinate system
%
%   UTMCoordinateSystem properties:
%   description - 'UTM' (Constant property)
%   zone - UTM zone
%
%   UTMCoordinateSystem methods:
%   xy - returns the UTM coordinates for the given lat lon coordinates
%
%   see also: ProjectedCoordinateSystem
    properties(Constant)
        description = 'UTM'
    end
    properties
        zone
    end
    methods
        function [x, y]=xy(obj,lat, lon)
        % returns the UTM Coordinates from lat lon coordinates
        %
        %   [x,y]=xy(obj, lat, lon) computes the x,y coordinates in the UTM
        %   projection given the lat,lon coordinates in the WGS84
        %   geographic coordinate system. If the zone property is empty, it
        %   is automatically determined and assigned. If the property is
        %   set, the given zone is used. 
        %
        %   see also: ProjectedCoordinateSystem, VMADCP
            if isempty(obj.zone)
                [x,y,obj.zone]=geo2utm(lat,lon);
            else
                [x,y]=geo2utm(lat,lon,obj.zone);
            end
        end
    end
    
end