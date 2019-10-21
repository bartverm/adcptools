classdef LatLonVisea < LatLonProvider
    methods
        function tf=has_data(obj,adcp)
            if isscalar(obj)
                tf=isfield(adcp,'VISEA_Extern') && all(isfield(adcp.VISEA_Extern,{'Latitudeseconds','Longitudeseconds'}));
            else
                tf=has_data@LatLonProvider(obj,adcp);
            end
        end
        function [lat, lon]=lat_lon(obj,adcp)
            if isscalar(obj)               
                if obj.has_data(adcp)
                   lat=reshape(adcp.VISEA_Extern.Latitudeseconds/3600,1,[]);
                   lon=reshape(adcp.VISEA_Extern.Longitudeseconds/3600,1,[]);
                else
                   lat=[];
                   lon=[];
                end
            else
                [lat, lon]=lat_lon@LatLonProvider(obj,adcp);
            end
        end
    end
end