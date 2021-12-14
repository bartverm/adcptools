classdef ADCPVerticalPositionFromWaterLevel < ADCPVerticalPosition
% Defines ADCP vertical position based on the water level
%
%   Subclasses should implement the get_water_level method
%
%   Constructor:
%   obj=ADCPVerticalPositionFromWaterLevel(...) constructs water level
%   object. Arguments passed to the constructor are handled according to
%   the type:
%   - WaterLevel objects are assigned to the water_level property
%   - double inputs are assigned to the depth_transducer property
%   - VMADCP objects are used to create a listener that updates the
%       water_level property, when the water_level_object property is
%       changed in the VMADCP object.
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
                elseif isa(arg,'VMADCP')
                    addlistener(arg,'water_level_object','PostSet',@obj.set_wl);
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
        function set_wl(obj, src, varargin)
            obj.water_level=src.water_level_object;
        end
    end
end