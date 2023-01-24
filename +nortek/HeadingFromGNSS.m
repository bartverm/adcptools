classdef HeadingFromGNSS < HeadingProvider
    methods (Access = protected)
        function val = get_has_data(~, adcp)
            val = isa(adcp, 'nortek.VMADCP');
        end
        function val = get_heading(obj, adcp)
            assert(isa(adcp,'nortek.VMADCP'),...
                'HeadingFromGNSS:NotVMADCP',...
                    'GNSS data only available for nortek.ADCP objects')
            val = interp1(adcp.gnss_time, adcp.gnss_heading, adcp.time,...
                'linear');
            if obj.heading_misalignment==0
                wrn_id = 'HeadingFromGNSS:ZeroMisalignment';
                warning(wrn_id,...
                    ['Heading misalignment set to zero. Make sure to ',...
                    'set it to the right value for proper heading. Not showing this warning again.'])
            end
        end
    end
end