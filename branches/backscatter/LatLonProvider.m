classdef LatLonProvider < matlab.mixin.Heterogeneous
    methods
        function tf=has_data(obj)
            tf = false;
            if numel(obj)>1
                for co=1:numel(obj)
                    tf=tf | obj(co).has_data;
                end
            end
        end
        function [lat, lon]=lat_lon(obj)
            lat=[];
            lon=[];
            co=1;
            while co <= numel(obj) && ~obj(co).has_data
                co=co+1;
            end
            if co > numel(obj)
                return
            else
                [lat,lon]=obj(co).lat_lon();
            end
        end
    end
end