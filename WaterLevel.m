classdef WaterLevel < handle
    methods (Abstract)
        wl=get_water_level(obj,time)
    end
end