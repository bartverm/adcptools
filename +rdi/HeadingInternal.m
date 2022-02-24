classdef HeadingInternal < HeadingProvider
    methods (Access=protected)
        function val=get_heading(~,adcp)
            val=double(adcp.raw.heading)/100;
        end
        function val=get_has_data(~,adcp)
            val = isfield(adcp.raw,'heading');
        end
    end
end