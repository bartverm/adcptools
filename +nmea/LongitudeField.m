classdef LongitudeField < nmea.LatLonField
    methods
        function obj = LongitudeField()
            obj = obj@nmea.LatLonField();
            obj.name = 'longitude';
            obj.pattern = obj.common_patterns.long;
            obj.format = ["%3f64" "%f64" "%c"];
        end
    end
end