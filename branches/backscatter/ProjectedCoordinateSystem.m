classdef ProjectedCoordinateSystem < handle
    properties (Constant, Abstract)
        description(1,:) char
    end
    methods (Abstract)
        [x,y]=xy(lat,lon)
    end
end