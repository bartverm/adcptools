classdef ShipVelocityFromBT < ShipVelocityProvider
    methods(Access=protected)
        function vel=get_ship_velocity(~,adcp,dst)
            vel=-adcp.btvel(dst);
        end
    end
end