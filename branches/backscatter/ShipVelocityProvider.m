classdef ShipVelocityProvider < matlab.mixin.Heterogeneous
    methods (Sealed)
        function vel=ship_velocity(obj,adcp,dst)
            validateattributes(adcp,{'VMADCP'},{'scalar'})
            validateattributes(dst,{'CoordinateSystem'},{'scalar'})
            if nargin < 3
                dst=CoordinateSystem.Earth;
            end
            if isscalar(obj)
                vel=obj.get_ship_velocity(adcp,dst);
            else
                vel=nan(1,adcp.nensembles,max(adcp.nbeams));
                for co=1:numel(obj)
                    isbad=isnan(vel);
                    cvel=obj(co).get_ship_velocity(adcp,dst);
                    vel(isbad)=cvel(isbad);
                end
            end
        end
    end
    methods (Abstract, Access=protected)
        vel=get_ship_velocity(obj,adcp,dst)
    end
end