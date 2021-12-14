classdef HeadingProviderInternal < HeadingProvider
% returns heading data from ADCP's internal compass
%
%   HeadingProviderInternal methods:
%   heading - returns heading
%   has_data - returns whether heading data are available
%
%   see also: HeadingProvider, HeadingProviderTFiles
    methods(Access=protected)
        function val=get_heading(~,adcp)
            val=double(adcp.raw.heading)/100;
        end
        function val=get_has_data(~,~)
            val=true;
        end
    end
end