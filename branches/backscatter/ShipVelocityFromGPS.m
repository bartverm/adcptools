classdef ShipVelocityFromGPS < ShipVelocityProvider
    methods(Access=protected)
        function vel=get_ship_velocity(~,adcp,dst)
            [x, y] = adcp.xy;
            t=adcp.time;
            dt=seconds(t(3:end)-t(1:end-2));
            dx=x(3:end)-x(1:end-2);
            dy=y(3:end)-y(1:end-2);
            vx=dx./dt;
            vy=dy./dt;
            vx=[vx(1) vx vx(end)];
            vy=[vy(1) vy vy(end)];
            vel=cat(3,vx,vy,zeros([size(vx),2]));
            vel=helpers.matmult(adcp.xform(dst, CoordinateSystem.Earth),vel);
        end
    end
end