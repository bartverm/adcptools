classdef ShipVelocityFromGPS < ShipVelocityProvider
% Computes ship velocity using GPS positions
%
%   Velocities are computed using GPS positions and the time the positions
%   where measured. Differentials of time and space are estimated using
%   central differences. Only horizontal ship velocity is computed.
%   Vertical ship velocity is set to zero.
%
%   ShipVelocityFromGPS methods:
%   ship_velocity - Ship velocity in m/s
%
%   see also: VMADCP, ProjectedCoordinateSystem, ShipVelocityProvider
    methods(Access=protected)
        function vel=get_ship_velocity(~,adcp,dst)
            [x, y] = deal(adcp.horizontal_position(1,:), adcp.horizontal_position(2,:));
            t=adcp.time;
            dt=seconds(t(3:end)-t(1:end-2));
            dx=x(3:end)-x(1:end-2);
            dy=y(3:end)-y(1:end-2);
            vx=dx./dt;
            vy=dy./dt;
            vx=[vx(1) vx vx(end)];
            vy=[vy(1) vy vy(end)];
            vel=cat(3,vx,vy,zeros([size(vx),2]));
            vel=helpers.matmult(adcp.xform(dst, CoordinateSystem.Earth,...
                'BottomTracking',true),vel);
            % last line transforms back to bottom_track system. This is to
            % match behaviour of systems with bottom tracking on different
            % beams than water tracking (e.g. sontek M9)
        end
    end
end