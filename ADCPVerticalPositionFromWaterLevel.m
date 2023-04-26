classdef ADCPVerticalPositionFromWaterLevel < ADCPVerticalPosition
% Defines ADCP vertical position based on the water level
%
%   Constructor:
%   obj=ADCPVerticalPositionFromWaterLevel(...)
%
%   WaterLevel methods:
%   get_water_level - returns the water level for each given time
%
%   see also: ADCPVerticalPosition
    properties
        depth_transducer (1,:) double {mustBeFinite, mustBeReal} = 0;
        water_level(1,1) WaterLevel = ConstantWaterLevel(0)
    end
    methods
        function val =  get.water_level(obj)
            warning(['water_level property is set for removal. Use the ',...
                'water_level_object property of the ADCP object to ',...
                'set the water level'])
            val = obj.water_level;
        end
        function set.water_level(obj, val)
            warning(['water_level property is set for removal. Use the ',...
                'water_level_object property of the ADCP object to ',...
                'set the water level'])
            obj.water_level = val;
        end
        function val =  get.depth_transducer(obj)
            warning(['depth_transducer property is set for removal. Use the ',...
                'depth_tranducer property of the VMADCP object to ',...
                'set the depth of the transducer'])
            val = obj.depth_transducer;
        end
        function set.depth_transducer(obj, val)
            warning(['depth_transducer property is set for removal. Use the ',...
                'depth transducer property of the VMADCP object to ',...
                'set the depth of the transducer'])
            obj.depth_transducer = val;
        end        
        function val=get_vertical_position(~,vmadcp)
%             if all(obj.depth_transducer==0)
%                 warning('Set depth of transducer for correct results')
%             end
            val=vmadcp.water_level - vmadcp.depth_transducer;
        end
    end
end