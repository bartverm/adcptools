classdef UTCField < nmea.Field
    methods
        function obj = UTCField()
            obj = obj@nmea.Field();
            obj.name = 'utc';
            obj.pattern = obj.common_patterns.time;
            obj.format = ["%2f32" "%2f32" "%f32"];
            obj.post_process = @(x) horzcat(x{:});
        end
    end
end