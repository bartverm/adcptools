classdef Bathymetry < handle
    properties
        water_level (1,1) WaterLevel=ConstantWaterLevel();
    end
    methods (Abstract)
        z=get_bed_elev(obj,x,y)
    end
    methods
        function d=get_depth(obj,x,y,time)
            if nargin > 3
                d=obj.get_bed_elev(x,y)+obj.water_level.get_water_level(time);
            else
                d=obj.get_bed_elev(x,y)+obj.water_level.get_water_level();
            end
        end
    end
end