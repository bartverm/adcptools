classdef HeadTiltFromAHRS < HeadTiltProvider
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            if ~isa(adcp, 'nortek.ADCP')
                error('HeadTiltFromAHRS only works with nortek.ADCP objects')
            end
            val = any(adcp.burst_has_data(adcp.BURST_BIT_AHRS));
        end
        function val = get_head_tilt_matrix(~, adcp)
            val = adcp.ahrs_matrix;
        end
    end
end