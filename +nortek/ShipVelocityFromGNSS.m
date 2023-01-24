classdef ShipVelocityFromGNSS < ShipVelocityProvider
    methods(Access=protected)
        function vel=get_ship_velocity(~,adcp,dst)
            if ~isa(adcp,'nortek.VMADCP')
                vel = nan(1,adcp.nensembles,3);
                return
            end
            t_gnss = adcp.gnss_time;
            t_adcp = adcp.time;
            u_gnss = adcp.gnss_vel_east;
            v_gnss = adcp.gnss_vel_north;
            w_gnss = adcp.gnss_vel_down;
            vel = cat(3,...
                interp1(t_gnss, u_gnss, t_adcp, 'linear'),...
                interp1(t_gnss, v_gnss, t_adcp, 'linear'),...
                interp1(t_gnss, w_gnss, t_adcp, 'linear')*0,... % for now not applying vertical velocity correction
                zeros(size(t_adcp)));
            vel=helpers.matmult(adcp.xform(dst, CoordinateSystem.Earth),vel);
        end
    end
end