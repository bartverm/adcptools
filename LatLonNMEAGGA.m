classdef LatLonNMEAGGA < LatLonProvider
% Reads geographic coordinates form GGA data in pd0 data structure
%
%   LatLonVisea methods:
%   has_data - whether the providers are able to give geographical data
%   lat_lon - get latitude and longitude
%
%   see also: VMADCP, LatLonNfilesGGA, LatLonProvider
    methods(Access=protected)
        function tf=get_has_data(~,adcp)
            tf= isfield(adcp.raw,'NMEAGGA') && ...
                all(isfield(adcp.raw.NMEAGGA,{'Lat','Long','deltaT'}));
        end
        function [lat, lon]=get_lat_lon(~,adcp)
           dt=adcp.raw.NMEAGGA.deltaT; % get time offset of GPS data
           fbad = ~isfinite(adcp.raw.NMEAGGA.Lat) |...
               ~isfinite(adcp.raw.NMEAGGA.Long) |...
               adcp.raw.NMEAGGA.Lat == 0 |...
               adcp.raw.NMEAGGA.Long == 0;
           dt(fbad) = nan;
           [~,fmin]=min(abs(dt),[],2,"omitnan"); % find closest GPS data for each ensemble
           fidx=sub2ind(size(adcp.raw.NMEAGGA.Lat),(1:size(adcp.raw.NMEAGGA.Lat,1))',fmin); % transform to linear index
           lat=reshape(adcp.raw.NMEAGGA.Lat(fidx),1,[]);
           lon=reshape(adcp.raw.NMEAGGA.Long(fidx),1,[]);
           fbad=lat==0 & lon==0;
           lat(fbad)=nan;
           lon(fbad)=nan;
        end
    end
end