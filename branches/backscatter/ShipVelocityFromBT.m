classdef ShipVelocityFromBT < ShipVelocityProvider
% Computes ship velocity from Bottom Tracking data
%
%   ShipVelocityFromBt methods:
%   ship_velocity - returns ship velocity in m/s
%
%   see also: VMADCP, ShipVelocityProvider, ShipVelocityFromGPS
    methods(Access=protected)
        function vel=get_ship_velocity(~,adcp,dst)
            vel=-adcp.btvel(dst);
        end
    end
end