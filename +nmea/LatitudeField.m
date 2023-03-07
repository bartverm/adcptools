classdef LatitudeField < nmea.LatLonField
    methods
        function obj = LatitudeField()
            obj = obj@nmea.LatLonField();
            obj.name = 'latitude';
            obj.pattern = obj.common_patterns.lat;
            obj.format = ["%2f64" "%f64" "%c"];
        end
    end
end