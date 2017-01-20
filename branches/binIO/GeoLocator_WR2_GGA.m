classdef (Abstract) GeoLocator_WR2_GGA < GeoLocator
    methods (Static)
        function [lat, long] = get_lat_long(adcp_structure)
            if ~isfield(adcp_structure,'gga_wr2_2')
                error('get_lat_long:no_wr2_2_gga','Could not find wr2 v.2 gga field in adcp structure')
            end
            nens=numel(adcp_structure.gga_wr2_2.dt);
            [lat, long]=deal(nan(nens,1));            
            for ce=1:nens
                ndat=numel(adcp_structure.gga_wr2_2.latitude{ce});
                if ndat < 2
                    continue
                end
                ew=ones(1,ndat);
                ew(adcp_structure.gga_wr2_2.lon_ew{ce}=='W')=-1;
                long(ce)=interp1(adcp_structure.gga_wr2_2.dt{ce},ew.*adcp_structure.gga_wr2_2.longitude{ce},0,'pchip');
                sn=ones(1,ndat);
                sn(adcp_structure.gga_wr2_2.lat_ns{ce}=='S')=-1;    
                lat(ce)=interp1(adcp_structure.gga_wr2_2.dt{ce},sn.*adcp_structure.gga_wr2_2.latitude{ce},0,'pchip');
            end
        end
    end
end