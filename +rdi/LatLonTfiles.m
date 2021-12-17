classdef LatLonTfiles < LatLonProvider
% Reads geographic coordinates form transect files
%
%   LatLonVisea methods:
%   has_data - whether the providers are able to give geographical data
%   lat_lon - get latitude and longitude
%
%   see also: VMADCP, LatLonNfilesGGA, LatLonProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf=isfield(adcp.raw,'tFiles') && ...
                all(isfield(adcp.raw.tFiles,{'lat','long'}));
        end
        function [lat, lon]=get_lat_lon(~,adcp)
           lat=reshape(adcp.raw.tFiles.lat,1,[]);
           lon=reshape(adcp.raw.tFiles.long,1,[]);
        end
    end
end