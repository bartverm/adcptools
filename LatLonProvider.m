classdef LatLonProvider < matlab.mixin.Heterogeneous
% Base class to obtain geographic coordinates from adcp data
%
%   This class can be subclassed to add provide a way to obtain geographic
%   coordinates from ADCP data. Different classes derived from this base
%   class can be combined in an array. When requesting the geographic
%   coordinates the class will query all the providers and the first that
%   is able to provide the data will be used.
%
%   LatLonProvider methods:
%   has_data - whether the providers are able to give geographical data
%   lat_lon - get latitude and longitude
%
%   see also: VMADCP, LatLonVisea, LatLonNfilesGGA

    methods (Sealed)
        function tf = has_data(obj, adcp)
        % returns whether geographic data is available
        %
        %   tf = has_data(obj, adcp) returns the array tf which has the
        %   same size as obj, with a logical value indicating whether the
        %   data providers are able to get latitude and longitude from the
        %   given adcp data. adcp must be an object of VMADCP class
        %
        %   see also: LatLonProvider, lat_lon
            validateattributes(adcp,{'VMADCP'},{'scalar'})
            tf=false(size(obj));
            for co=1:numel(obj)
                tf(co)= obj(co).get_has_data(adcp);
            end
        end
        function [lat, lon] = lat_lon(obj, adcp)
        % return latitude and longitude from adcp data
        %
        %   [lat, lon] = lat_lon(obj, adcp) returns latitude and longitude
        %   information using the first provider in obj that can provide
        %   these data from the adcp object of class VMADCP.
        %
        %   see also: LatLonProvider, has_data
            validateattributes(adcp,{'VMADCP'},{'scalar'})
            fprovider = find(obj.has_data(adcp), 1);
            if isempty(fprovider)
                [lat, lon]=deal(nan(1, adcp.nensembles));
            else
                [lat,lon]=obj(fprovider).get_lat_lon(adcp);
            end
        end
    end
    methods (Access=protected, Abstract) % methods to implement in subclasses
        get_lat_lon(adcp) 
        get_has_data(adcp)
    end
end