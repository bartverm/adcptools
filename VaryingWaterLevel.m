classdef VaryingWaterLevel < WaterLevel
    properties
        water_level (1,:) double {mustBeReal, mustBeFinite} = 0
    end
    methods
        function obj=VaryingWaterLevel(varargin)
            for ca=1:nargin
                if isa(varargin{ca},'double')
                    obj.water_level=varargin{ca};
                else
                    warning(['Unused input number ',num2str(ca),' of type ',class(varargin{ca})])
                end
            end
        end
        function val=get_water_level(obj,adcp)
            assert(numel(obj.water_level)==adcp.nensembles,'Elements in water level series do not match number of ensembles')
            val=obj.water_level;
        end
        
    end
end