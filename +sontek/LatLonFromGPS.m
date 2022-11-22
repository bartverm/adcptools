classdef LatLonFromGPS < LatLonProvider
    methods(Access = protected)
        function [lat, lon] = get_lat_lon(~, adcp) 
            lat = adcp.raw.GPS.Latitude';
            lon = adcp.raw.GPS.Longitude';
        end
        function val = get_has_data(~, adcp)
            val = isfield(adcp.raw,'GPS') &&...
                isfield(adcp.raw.GPS,'Longitude') &&...
                isfield(adcp.raw.GPS,'Latitude');
        end
    end
end