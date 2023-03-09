classdef WaterLevel < handle
% Defines ADCP vertical position based on the water level
%
%   Subclasses should implement the get_water_level method
%
%   WaterLevel methods:
%   get_water_level - returns the water level for each given time
%   get_depth - compute the depth for a given time and elevation
%
%   see also: ConstantWaterLevel, VaryingWaterLevel
    methods (Abstract)
        wl=get_water_level(obj,time)
    end
    methods
        function depth=get_depth(obj,z,time)
            % Compute the depth for a given time and elevation
            %
            %   depth = get_depth(obj, z, time) get the depth for the elevation z at
            %       the given time
            %
            %   see also: VaryingWaterLevel
            wl=obj.get_water_level(time);
            depth=wl-z;
        end
    end
end