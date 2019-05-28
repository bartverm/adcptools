classdef WaterLevel < handle
    methods (Abstract)
        wl=get_water_level(obj,time)
    end
    methods
        function depth=get_depth(obj,z,time)
            if nargin > 2
                wl=obj.get_water_level(time);
            else
                wl=obj.get_water_level();
            end
            depth=wl-z;
        end
    end
end