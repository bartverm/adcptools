classdef TiltsInternal < TiltsProvider
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            val = adcp.has_tilts_internal;
        end
        function val = get_pitch(~, adcp)
            val = adcp.pitch_internal;
        end
        function val = get_roll(~, adcp)
            val = adcp.roll_internal;
        end
    end
end