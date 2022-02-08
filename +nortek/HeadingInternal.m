classdef HeadingInternal < HeadingProvider
    methods (Access = protected)
        function val = get_has_data(~, adcp)
            val = adcp.has_heading_internal;
        end
        function val = get_heading(~, adcp)
            val = adcp.heading_internal;
        end
    end
end