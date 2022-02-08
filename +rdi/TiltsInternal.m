classdef TiltsInternal < TiltsProvider
    methods (Access=protected)
        function pitch=get_pitch(~,adcp)
            pitch=double(adcp.raw.pitch)/100;
            pitch=atand(tand(pitch).*cosd(adcp.roll));
        end
        function roll=get_roll(~,adcp)
            roll=double(adcp.raw.roll)/100;
        end
        function val=get_has_data(~,adcp)
            val = isfield(adcp.raw,'pitch') & isfield(adcp.raw,'roll');
            val = repmat(val, 1, adcp.nensembles);
        end
    end
end