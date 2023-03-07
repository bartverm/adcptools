classdef RMCMessage < nmea.Message
    methods
        function obj = RMCMessage()
            obj.msg_id_pattern = 'RMC';
            obj.name = obj.msg_id_pattern;
            cp = nmea.Field.common_patterns;
            obj.fields = [
                nmea.UTCField;
                nmea.Field('status',"%1c",cp.status);
                nmea.LatitudeField;
                nmea.LongitudeField;
                nmea.Field('ground_speed',"%f32", cp.float);
                nmea.Field('track_made_good',"%f32", cp.float);
                nmea.Field('date', ["%2f32", "%2f32", "%2f32"],...
                    cp.ndgt(6),...
                    @(x) [x{:}]);
                nmea.Field('magnetic_variation', ["%f32", "%1s"],...
                    [cp.float,',[EW]?'],...
                    @(x) (strcmp(x{2},'W')*2-1).*x{1});
                nmea.Field('mode',"%1c",'[ADEMSN]?');
            ];
        end
    end
end