classdef HeadingFromAHRS < HeadingProvider
    methods (Access = protected)
        function val = get_has_data(~, adcp)
            if ~isa(adcp, 'nortek.ADCP')
                error('HeadingFromAHRS only works with nortek.ADCP objects')
            end
            val = any(adcp.burst_has_data(nortek.BurstBit.AHRS));
        end
        function val = get_heading(~, adcp)
            tm_ahrs = adcp.ahrs_matrix;
            val = 90+atan2d(-tm_ahrs(:,:,2,1), tm_ahrs(:,:,1,1));
            val (val < 0 ) = val (val < 0 ) +360;
            val (val >= 360 ) = val (val >=360 ) - 360;
        end
    end
end