classdef HeadingProviderHFiles < HeadingProvider
% returns heading form h-files
%
%   HeadingProviderHFiles methods:
%   heading - returns heading
%   has_data - returns whether heading data are available
%
%   see also: HeadingProviderInternal, HeadingProvider
    methods(Access=protected)
        function val=get_heading(~,adcp)
            val=reshape(double(adcp.raw.hFiles.HDT.heading),1,[]);
        end
        function val=get_has_data(~,adcp)
            val=isfield(adcp.raw,'hFiles') &&...
                isfield(adcp.raw.hFiles,'HDT') && ...
                isfield(adcp.raw.hFiles.HDT,'heading');
        end
    end
end