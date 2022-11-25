classdef TiltsInternal < TiltsProvider
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            val = isfield(adcp.raw,'Compass') &&...
                isfield(adcp.raw.Compass, 'Pitch') &&...
                isfield(adcp.raw.Compass, 'Roll') && ...
                isfield(adcp.raw.Compass, 'Units') && ...
                isfield(adcp.raw.Compass.Units, 'Pitch') &&...
                isfield(adcp.raw.Compass.Units, 'Roll');
        end
        function val = get_pitch(~, adcp)
            val = sontek.TiltsInternal.get_tilt(adcp,'Pitch');
        end
        function val = get_roll(~, adcp)
            val = sontek.TiltsInternal.get_tilt(adcp,'Roll');
        end
    end
    methods(Access=protected, Static)
        function val = get_tilt(adcp,name)
            val = adcp.raw.Compass.(name)(:,1)';
            unit = adcp.raw.Compass.Units.(name)(adcp.raw.file_id);
            is_rad = strcmp(unit,'rad');
            val(is_rad) = rad2deg(val(is_rad));
            is_deg = strcmp(unit,'deg');
            is_known = is_deg | is_rad;
            assert(all(is_known),...
                ['Not all ', name, ' angle units are known'])
        end
    end
end