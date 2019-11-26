classdef LatLonGGA < LatLonProvider
% Reads geographic coordinates form GGA data in pd0 data structure
%
%   LatLonVisea methods:
%   has_data - whether the providers are able to give geographical data
%   lat_lon - get latitude and longitude
%
%   see also: VMADCP, LatLonNfilesGGA, LatLonProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf= isfield(adcp.raw,'GGA') && ...
                all(isfield(adcp.raw.GGA,{'lat','long'}));
        end
        function [lat, lon]=get_lat_lon(~,adcp)
           lat=reshape(adcp.raw.GGA.lat,1,[]);
           lon=reshape(adcp.raw.GGA.long,1,[]);
        end
    end
end