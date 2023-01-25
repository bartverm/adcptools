classdef LatLonFromGNSS < LatLonProvider
    methods (Access=protected)
        function [lat, lon] = get_lat_lon(~, adcp)
            lat = adcp.gnss_latitude;
            lon = adcp.gnss_longitude;
            gps_time = adcp.gnss_time;
            adcp_time = adcp.time;
            lat = interp1(gps_time, lat, adcp_time,'nearest');
            lon = interp1(gps_time, lon, adcp_time,'nearest');
        end
        function val = get_has_data(~, adcp)
            val = isprop(adcp,'gnss_latitude') &&...
                isprop(adcp,'gnss_longitude');
        end
    end
end