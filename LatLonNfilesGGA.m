classdef LatLonNfilesGGA < LatLonProvider
% Reads geographic coordinates form GGA string in n-files
%
%   LatLonVisea methods:
%   has_data - whether the providers are able to give geographical data
%   lat_lon - get latitude and longitude
%
%   see also: VMADCP, LatLonNfilesGGA, LatLonProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf=isfield(adcp.raw,'nFiles') && ...
                isfield(adcp.raw.nFiles,'GGA') && ...
                all(isfield(adcp.raw.nFiles.GGA,{'lat','long'}));
        end
        function [lat, lon]=get_lat_lon(~,adcp)
           lat=reshape(adcp.raw.nFiles.GGA.lat,1,[]);
           lon=reshape(adcp.raw.nFiles.GGA.long,1,[]);
        end
    end
end