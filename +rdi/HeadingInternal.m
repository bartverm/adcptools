classdef HeadingInternal < HeadingProvider
    methods (Access=protected)
        function val=get_heading(~,adcp)
            val=double(adcp.raw.heading)/100;
        end
        function val=get_has_data(~,adcp)
            val=all(strcmp(adcp.raw.sensource(:,5),'1'));
        end
    end
end