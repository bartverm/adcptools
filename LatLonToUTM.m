classdef LatLonToUTM < LatLonToProjection
% UTM coordinate system
%
%   LatLonToUTM properties:
%   description - 'UTM' (Constant property)
%   zone - UTM zone
%
%   LatLonToUTM methods:
%   xy - returns the UTM coordinates for the given lat lon coordinates
%
%   see also: LatLonToProjection
    properties
        description = 'UTM'
    end
    properties
        zone
    end
    methods
        function obj=LatLonToUTM(varargin)
            obj=obj@LatLonToProjection(varargin{:});
        end
        function [x, y]=xy(obj,lat, lon)
        % returns the UTM Coordinates from lat lon coordinates
        %
        %   [x,y]=xy(obj, lat, lon) computes the x,y coordinates in the UTM
        %   projection given the lat,lon coordinates in the WGS84
        %   geographic coordinate system. If the zone property is empty, it
        %   is automatically determined and assigned. If the property is
        %   set, the given zone is used. 
        %
        %   see also: LatLonToProjection, VMADCP
            if isempty(obj.zone)
                [x,y,obj.zone]=helpers.geo2utm(lat,lon);
            else
                [x,y]=helpers.geo2utm(lat,lon,obj.zone);
            end
        end
    end
    
end