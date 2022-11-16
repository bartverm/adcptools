classdef ErrorVelocityFilter < Filter
   properties
       max_err_vel(1,1) double {mustBeNonnegative} = 0.8;
   end
   methods
        function obj=ErrorVelocityFilter()
            obj.description='Error velocity filter';
        end
    end
   methods (Access=protected)
       function bad=bad_int(obj,adcp)
           if isscalar(obj)
               dst=adcp.coordinate_system;
               dst(dst<2)=CoordinateSystem.Instrument;
               vel=adcp.velocity(dst, Filter); % Get unfiltered velocity with in at least instrument coordinates, to get error velocity
               bad=repmat(vel(:,:,4)>=obj.max_err_vel,[1 1 4]);
           else
               bad=bad_int@Filter(obj,adcp);
           end
       end
   end
end