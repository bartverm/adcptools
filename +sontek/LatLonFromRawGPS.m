classdef LatLonFromRawGPS < LatLonProvider
    methods(Access = protected)
        function [lat, lon] = get_lat_lon(~, adcp) 
            lat = adcp.raw.RawGPSData.GgaLatitude';
            lon = adcp.raw.RawGPSData.GgaLongitude';
            lat(lat==0) = nan;
            lon(lon==0) = nan;
            lat = mean(lat,1,'omitnan');
            lon = mean(lon,1,'omitnan');
        end
        function val = get_has_data(~, adcp)
            val = isfield(adcp.raw,'RawGPSData') &&...
                isfield(adcp.raw.RawGPSData,'GgaLongitude') &&...
                isfield(adcp.raw.RawGPSData,'GgaLatitude') &&...
                any(adcp.raw.RawGPSData.GgaLatitude~=0,"all") && ...
                any(adcp.raw.RawGPSData.GgaLongitude~=0,"all");
        end
    end
end