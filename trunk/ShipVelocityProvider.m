classdef ShipVelocityProvider < matlab.mixin.Heterogeneous
% Base class to implement methods that yield the ship's velocity
%
%   ShipVelocityProvider methods:
%   ship_velocity - velocity of the ship in m/s
%
%   see also: VMADCP, ShipVelocityFromBT, ShipVelocityFromGPS

    methods (Sealed)
        function vel=ship_velocity(obj,adcp,dst)
% Computes ship velocity in m/s
%
%   vel=ship_velocity(obj, adcp, dst) computes the velocity of the ship in
%   m/s. This velocity is used to obtain water velocity from the apparent
%   velocity measured by the ADCP, for a moving ADCP.
%   Inputs:
%   obj: vector of ShipVelocityProvider objects
%   adcp: scalar VMADCP object
%   dst: CoordinateSystem object defining destination coordinate system
%
%   When multiple ShipVelocityProvider objects are passed, the function
%   will use the first provider to obtain the velocity and then proceed to
%   the next provider to fill the missing values, until the last provider
%   given.
%  
%   see also: VMADCP, CoordinateSystem
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
        vel=get_ship_velocity(obj, adcp, dst)
        % Abstract method computing the ship velocity from a scalar
        % ShipVelocityProvider subclass object. This method must be
        % implemented in a subclass
    end
end