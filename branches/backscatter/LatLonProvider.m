classdef LatLonProvider < matlab.mixin.Heterogeneous
    methods (Sealed)
        function tf=has_data(obj,adcp)
            tf=false(size(obj));
            for co=1:numel(obj)
                tf(co)= obj(co).get_has_data(adcp);
            end
        end
        function [lat, lon]=lat_lon(obj,adcp)
            fprovider=find(obj.has_data(adcp),1);
            if isempty(fprovider)
                [lat, lon]=deal(nan(1,adcp.nensembles));
            else
                [lat,lon]=obj(fprovider).get_lat_lon(adcp);
            end
        end
    end
    methods (Access=protected, Abstract)
        get_lat_lon(adcp)
        get_has_data(adcp)
    end
end