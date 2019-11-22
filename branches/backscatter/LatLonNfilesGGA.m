classdef LatLonNfilesGGA < LatLonProvider
    methods(Access=protected)
        function tf=get_has_data(obj,adcp)
            tf=isfield(adcp,'nFiles') && ...
                isfield(adcp.nFiles,'GGA') && ...
                all(isfield(adcp.nFiles.GGA,{'lat','long'}));
        end
        function [lat, lon]=get_lat_lon(obj,adcp)
           lat=reshape(adcp.nFiles.GGA.lat,1,[]);
           lon=reshape(adcp.nFiles.GGA.long,1,[]);
        end
    end
end