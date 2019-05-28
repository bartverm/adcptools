classdef ConstantWaterLevel < WaterLevel
    properties
        level=0;
    end
    methods
        function obj=ConstantWaterLevel(varargin)
            if nargin > 0
                obj.level=varargin{1};
            end
        end
        function set.level(obj,in_level)
            validateattributes(in_level,{'numeric'},{'scalar','real'})
            obj.level=in_level;
        end
        function wl=get_water_level(obj,time)
            if nargin>1
                wl=ones(size(time))*obj.level;
            else
                wl=obj.level;
            end
        end
    end
end