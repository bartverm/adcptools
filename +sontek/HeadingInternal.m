classdef HeadingInternal < HeadingProvider
    methods(Access = protected)
        function val = get_has_data(~, adcp)
            val = isfield(adcp.raw,'System') &&...
                isfield(adcp.raw.System,'True_North_ADP_Heading') &&...
                isfield(adcp.raw.System,'Units') &&...
                isfield(adcp.raw.System.Units,'Heading');
        end
        function val = get_heading(~, adcp)
            val = reshape(adcp.raw.System.True_North_ADP_Heading,1,[]);
            unit = adcp.raw.System.Units.Heading(adcp.raw.file_id);
            is_deg = strcmp(unit,'deg');
            is_rad = strcmp(unit,'rad');
            val(is_rad) = rad2deg(val(is_rad));
            val = val - 90;
            is_known = is_deg | is_rad;
            assert(all(is_known),'Not all heading angle units are known')
            if isfield(adcp.raw,'Setup') &&...
                    isfield(adcp.raw.Setup,'magneticDeclination')
                val = val +...
                    reshape(adcp.raw.Setup.magneticDeclination(...
                    adcp.raw.file_id),1,[]);
            end
        end
    end
end