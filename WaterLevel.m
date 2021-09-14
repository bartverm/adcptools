classdef WaterLevel < ADCPVerticalPosition
% Abstract ADCPVerticalPosition class defining the water level
%
%   Subclasses should implement the get_water_level method
%
%   WaterLevel methods:
%   get_water_level - returns the water level for each given time
%   get_depth - compute the depth for a given time and elevation
%
%   see also: ConstantWaterLevel
    properties
        depth_transducer (1,:) double {mustBeFinite, mustBeReal} = 0;
    end
    methods (Abstract)
        wl=get_water_level(obj,adcp)
    end
    methods
        function val=get_vertical_position(obj,adcp)
            if all(obj.depth_transducer==0)
                warning('Set depth of transducer for correct results')
            end
            val=obj.get_water_level(adcp)-obj.depth_transducer;
        end
        function depth=get_depth(obj,z,adcp)
% Compute the depth for a given time and elevation
%
%   depth = get_depth(obj, z, time) get the depth for the elevation z at
%       the given time
%
%   see also: WaterLevel
            if nargin > 2
                wl=obj.get_water_level(adcp);
            else
                wl=obj.get_water_level();
            end
            depth=wl-z;
        end
    end
end