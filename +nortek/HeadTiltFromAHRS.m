classdef HeadTiltFromAHRS < nortek.HeadTiltProvider
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            if ~isa(adcp, 'nortek.ADCP')
                error('HeadTiltFromAHRS only works with nortek.ADCP objects')
            end
            val = any(adcp.burst_has_data(nortek.BurstBit.AHRS));
        end
        function val = get_head_tilt_matrix(~, adcp)
            val = adcp.ahrs_matrix;
            val = cat(4,...
                cat(3,...
                val,...
                zeros([size(val,[1,2]),1,3])),...
                cat(3,...
                zeros([size(val,[1,2]),3,1]),...
                ones([size(val,[1,2]),1,1])));
        end
    end
end