classdef ADCPFixedHorizontalPosition < ADCPHorizontalPosition
% Defines the ADCP position as a fixed position (moored deployment)
    properties
        description = '';
        position (2,1) double {mustBeFinite, mustBeReal} = [0; 0];
    end
    methods
        function obj=ADCPFixedPosition(pos)
            if nargin > 0
                obj.position=pos;
            end
        end
    end
    methods(Access=protected)
        function tf=get_has_data(~,~)
            tf=true;
        end
        function pos=get_horizontal_position(obj,adcp)
            validateattributes(adcp,{'ADCP'},{'scalar'});
            pos=repmat(obj.position,[1,adcp.nensembles]);
        end
    end
end