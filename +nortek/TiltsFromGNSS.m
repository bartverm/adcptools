classdef TiltsFromGNSS < TiltsProvider
    properties
        heading_misalignment = 0;
        roll_misalignment = 0;
        pitch_misalignment = 0;
    end
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            val = isa(adcp, 'nortek.VMADCP');
        end
        function val = get_pitch(obj, adcp)
            obj.check_vmadcp(adcp)
            val = interp1(adcp.gnss_time, adcp.gnss_pitch, adcp.time,...
                'linear');
            val = obj.correct_bounds(val + obj.pitch_misalignment);
        end
        function val = get_roll(obj, adcp)
            obj.check_vmadcp(adcp)
            val = interp1(adcp.gnss_time, adcp.gnss_roll, adcp.time,...
                'linear');
            val = obj.correct_bounds(val + 180 + obj.roll_misalignment);
        end
        function check_vmadcp(~,adcp)
            assert(isa(adcp,'nortek.VMADCP'),...
                'HeadingFromGNSS:NotVMADCP',...
                'GNSS data only available for nortek.ADCP objects')
        end
    end
    methods
        function set_tilt_misalignment(obj, adcp)
            obj.roll_misalignment = 0;
            obj.pitch_misalignment = 0;
            AHRS = nortek.TiltsFromAHRS;
            roll_ahrs = AHRS.roll(adcp);
            pitch_ahrs = AHRS.pitch(adcp);
            roll_gnss = obj.roll(adcp);
            pitch_gnss = obj.pitch(adcp);
            obj.roll_misalignment = obj.dangle(roll_ahrs, roll_gnss);
            obj.pitch_misalignment = obj.dangle(pitch_ahrs, pitch_gnss);
        end
    end
    methods(Static)
        function val = correct_bounds(val)
            val(val > 180) = val(val > 180) - 360;
            val(val < -180) = 360 + val(val < -180);
        end
        function val = dangle(val1, val2)
            val = val1 - val2;
            val = mod(val + 180, 360) - 180;
            val = mean(val,'all','omitnan');
        end
    end
end
