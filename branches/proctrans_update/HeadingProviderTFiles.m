classdef HeadingProviderTFiles < HeadingProvider
    methods(Access=protected)
        function val=get_heading(~,adcp)
            val=adcp.raw.tFiles.corrHead;
        end
        function val=get_has_data(~,adcp)
            val=isfield(adcp.raw,'tFiles') && isfield(adcp.raw.tFiles,'corrHead');
        end
    end
end