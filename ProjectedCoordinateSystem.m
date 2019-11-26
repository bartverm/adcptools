classdef ProjectedCoordinateSystem < handle
% Base class for geographic to projected coordinate system conversion
%
%   ProjectedCoordinateSystem properties
%   description - holds an abbreviation for the coordinate system
%
%   ProjectedCoordinateSystem methods
%   xy - returns projected coordinates for the given lat,lon coordinates
%
%   see also: UTMCoordinateSystem, RDCoordinateSystem
    properties(Constant, Abstract)
        description
    end
    methods (Abstract)
        [x,y]=xy(lat,lon)
    end
end