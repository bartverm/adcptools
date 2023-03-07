classdef LatLonField < nmea.Field
    methods
        function obj = LatLonField()
            obj.post_process =...
                @(x) ((x{3}=='N' | x{3}=='E')*2-1).*(x{1}+x{2}/60);
        end
    end
end