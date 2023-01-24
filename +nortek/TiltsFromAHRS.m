classdef TiltsFromAHRS < TiltsProvider
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            if ~isa(adcp, 'nortek.ADCP')
                error('TiltFromAHRS only works with nortek.ADCP objects')
            end
            val = any(adcp.burst_has_data(nortek.BurstBit.AHRS));
        end
        function val = get_pitch(obj, adcp)
            tm_ahrs = adcp.ahrs_matrix;
            ahrs_roll = obj.get_roll(adcp);
            val = atan2d(tm_ahrs(:,:,3,1), tm_ahrs(:,:,3,2)./sind(ahrs_roll));
        end
        function val = get_roll(~, adcp)
            tm_ahrs = adcp.ahrs_matrix;
            val = atan2d(tm_ahrs(:,:,3,2), tm_ahrs(:,:,3,3));
        end
    end
end
