classdef SkipField < nmea.Field
    methods
        function obj=SkipField(pattern)
            obj = obj@nmea.Field();
            obj.name = 'skip';
            obj.format = "";
            if nargin > 0
                obj.pattern = pattern;
            end
        end
    end
    methods(Access = protected)
        function val = get_named_pattern(obj)
            val = obj.pattern + "?";
        end
    end
end