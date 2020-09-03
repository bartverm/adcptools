classdef Bathymetry < handle
    properties
        water_level (1,1) WaterLevel=ConstantWaterLevel;
    end
    methods (Abstract)
        z=get_bed_elev(obj,x,y)
    end
    methods
        function obj=Bathymetry(varargin)
            for count_var=1:nargin
                cvar=varargin{count_var};
                if isa(cvar,'WaterLevel')
                    obj.water_level=cvar;
                end
            end
        end
        function d=get_depth(obj,x,y,time)
            if nargin > 3
                d=obj.water_level.get_depth(obj.get_bed_elev(x,y),time);
            else
                d=obj.obj.water_level.get_depth(obj.get_bed_elev(x,y));
            end
        end
        function plot(obj,varargin)
            obj(1).plot;
            hold_stat=get(gca,'NextPlot');
            hold on;
            if ~isscalar(obj)
                for ce=2:numel(obj)
                    obj(ce).plot(varargin{:})
                end
            end
            set(gca,'NextPlot',hold_stat);
        end
    end
end