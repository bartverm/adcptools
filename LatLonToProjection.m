classdef LatLonToProjection < ADCPHorizontalPosition
% Base class for geographic to projected coordinate system conversion
%
%   LatLonToProjection properties
%   description - holds an abbreviation for the coordinate system
%
%   LatLonToProjection methods
%   xy - returns projected coordinates for the given lat,lon coordinates
%
%   see also: UTMCoordinateSystem, RDCoordinateSystem
    properties
    % VMADCP/ll_provider 
    %
    %   Latitude, Longitude provider. Array of LatLonProvider objects. The
    %   first provider that has valid LatLon data is used. Default is
    %   [LatLonVisea; LatLonNfiles; LatLonTfiles; LatLonNMEAGGA; LatLonGGA]
    %
    %   see also: VMADCP, LatLonProvider
        ll_provider (:,1) LatLonProvider    
    end
    methods
        function obj=LatLonToProjection(varargin)
            for ca = 1 : nargin
                if isa(varargin{ca},'LatLonProvider')
                    obj.ll_provider = varargin{ca};
                end
            end
        end
    end
    methods(Access=protected)
        function tf=get_has_data(obj,adcp)
            tf=any(obj.ll_provider.has_data(adcp));
        end
        function pos=get_horizontal_position(obj,adcp)
            [lat,lon]=obj.ll_provider.lat_lon(adcp);
            [x,y]=obj.xy(lat,lon);
            pos=[x;y];
        end
    end
    methods (Abstract)
        [x,y]=xy(lat,lon)
    end
end