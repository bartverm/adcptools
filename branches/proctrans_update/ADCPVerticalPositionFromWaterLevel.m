classdef ADCPVerticalPositionFromWaterLevel < ADCPVerticalPosition
% Defines ADCP vertical position based on the water level
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
        water_level(1,1) WaterLevel = ConstantWaterLevel(0)
    end
    methods
        function obj=ADCPVerticalPositionFromWaterLevel(varargin)
            for ca=1:numel(varargin)
                arg=varargin{ca};
                if isa(arg, 'WaterLevel')
                    obj.water_level=arg;
                elseif isa(arg,'double')
                    obj.depth_transducer=arg;
                else
                    warning(['Unhandled input of type: ', class(arg)])
                end
            end
        end
        function val=get_vertical_position(obj,adcp)
%             if all(obj.depth_transducer==0)
%                 warning('Set depth of transducer for correct results')
%             end
            val=obj.water_level.get_water_level(adcp.time)-obj.depth_transducer;
        end
    end
end