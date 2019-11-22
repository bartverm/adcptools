classdef LatLonVisea < LatLonProvider
    methods(Access=protected)
        function tf=get_has_data(obj,adcp)
            tf=isfield(adcp,'VISEA_Extern') && all(isfield(adcp.VISEA_Extern,{'Latitudeseconds','Longitudeseconds'}));
        end
        function [lat, lon]=get_lat_lon(obj,adcp)
            lat=reshape(adcp.VISEA_Extern.Latitudeseconds/3600,1,[]);
            lon=reshape(adcp.VISEA_Extern.Longitudeseconds/3600,1,[]);
        end
    end
end