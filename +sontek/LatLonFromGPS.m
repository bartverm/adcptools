classdef LatLonFromGPS < LatLonProvider
    methods(Access = protected)
        function [lat, lon] = get_lat_lon(~, adcp) 
            lat = adcp.raw.GPS.Latitude';
            lon = adcp.raw.GPS.Longitude';
            lat(lat==0) = nan;
            lon(lon==0) = nan;
        end
        function val = get_has_data(~, adcp)
            val = isfield(adcp.raw,'GPS') &&...
                isfield(adcp.raw.GPS,'Longitude') &&...
                isfield(adcp.raw.GPS,'Latitude') &&...
                any(adcp.raw.GPS.Latitude~=0) && ...
                any(adcp.raw.GPS.Longitude~=0);
        end
    end
end