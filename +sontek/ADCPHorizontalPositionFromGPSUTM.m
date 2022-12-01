classdef ADCPHorizontalPositionFromGPSUTM < ADCPHorizontalPosition
    properties
        description = 'UTM';
    end
    methods(Access = protected)
        function pos = get_horizontal_position(~, adcp) 
            pos = adcp.raw.GPS.UTM';
            pos(pos==0) = nan;
        end
        function val = get_has_data(~, adcp)
            val = isfield(adcp.raw,'GPS') &&...
                isfield(adcp.raw.GPS,'UTM') &&...
                any(adcp.raw.GPS.UTM~=0,"all");
        end
    end
end