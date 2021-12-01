classdef HeadingProviderInternal < HeadingProvider
    methods(Access=protected)
        function val=get_heading(~,adcp)
            val=double(adcp.raw.heading)/100;
        end
        function val=get_has_data(~,~)
            val=true;
        end
    end
end