classdef LatLonVisea < LatLonProvider
% Reads geographic coordinates form VISEA extern files
%
%   LatLonVisea methods:
%   has_data - whether the providers are able to give geographical data
%   lat_lon - get latitude and longitude
%
%   see also: VMADCP, LatLonNfilesGGA, LatLonProvider

    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf=isfield(adcp.raw,'VISEA_Extern') && all(isfield(adcp.raw.VISEA_Extern,{'Latitudeseconds','Longitudeseconds'}));
        end
        function [lat, lon]=get_lat_lon(~,adcp)
            lat=reshape(adcp.raw.VISEA_Extern.Latitudeseconds/3600,1,[]);
            lon=reshape(adcp.raw.VISEA_Extern.Longitudeseconds/3600,1,[]);
        end
    end
end