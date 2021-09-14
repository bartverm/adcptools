classdef ADCPFixedVerticalPosition < ADCPVerticalPosition
% Defines the ADCP position as a fixed position (moored deployment)
    properties
        position (1,1) double {mustBeFinite, mustBeReal} = 0;
    end
    methods
        function obj=ADCPFixedPosition(pos)
            if nargin > 0
                obj.position=pos;
            end
        end
        function pos=get_vertical_position(obj,adcp)
            validateattributes(adcp,{'ADCP'},{'scalar'});
            pos=ones(1,adcp.nensembles)*obj.position;
        end
    end
end