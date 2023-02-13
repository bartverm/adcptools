classdef ADCPHorizontalPosition < matlab.mixin.Heterogeneous & handle
    % Abstract class defining the horizontal position of the ADCP
    %
    %   ADCPHorizontalPosition methods:
    %   get_horizontal_position - returns the ADCP position for the given time
    %
    % see also: ADCPFixedHorizontalPosition
    properties(Abstract)
        description
    end
    methods(Sealed)
        function tf = has_data(obj, adcp)
% returns whether geographic data is available
%
%   tf = has_data(obj, adcp) returns the array tf which has the
%   same size as obj, with a logical value indicating whether the
%   data providers are able to get latitude and longitude from the
%   given adcp data. adcp must be an object of VMADCP class
%
%   see also: LatLonProvider, lat_lon
            validateattributes(adcp,{'ADCP'},{'scalar'})
            tf=false(size(obj));
            for co=1:numel(obj)
                tf(co)= obj(co).get_has_data(adcp);
            end
        end
        function pos = horizontal_position(obj,adcp)
        % ADCPPosition/horizontal_position
        %
        % Method should be implemented by subclasses and return a 2xN
        % matrix given an 1XN datetime row vector and adcp object. The rows
        % of the output define the x,y position of the ADCP in m
        %
        % see also: ADCPHorizontalPosition
            fprovider = find(obj.has_data(adcp), 1);
            if isempty(fprovider)
                pos=deal(nan(2, adcp.nensembles));
            else
                pos=obj(fprovider).get_horizontal_position(adcp);
            end
        end
    end
    methods (Abstract, Access=protected)

        pos=get_horizontal_position(obj, adcp)
        tf=get_has_data(obj,adcp)
    end
end