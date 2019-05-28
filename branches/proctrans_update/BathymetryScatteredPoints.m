classdef BathymetryScatteredPoints < Bathymetry
    properties
        pos (3,:) double {mustBeFinite} = double.empty(3,0);
        interpolator (1,1) Interpolator = LoessInterpolator();
    end
    methods
        function z=get_bed_elev(obj,qpos)
            validateattributes(qpos,{'double'},{'2d','finite'})
            obj.interpolator.known=obj.pos;
            z=obj.interpolator.interpolate(qpos);
        end
        function plot3(obj)
           scatter3(reshape(obj.pos(1,:,:),1,[]),reshape(obj.pos(2,:,:),1,[]),reshape(obj.pos(3,:,:),1,[]),13,reshape(obj.pos(3,:,:),1,[]),'filled') 
        end
    end
end