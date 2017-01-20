classdef (Abstract) GeoLocator < handle
    methods (Abstract, Static)
        [lat, long] = get_lat_long(adcp_structure);
    end
end