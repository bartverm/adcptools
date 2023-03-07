classdef RDENSMessage < nmea.Message
    methods
        function obj=RDENSMessage()
            obj.talker_id_pattern = 'RD';
            obj.msg_id_pattern = 'ENS';
            obj.name = 'RDENS';
            obj.needs_valid_checksum = false;
            obj.fields = [nmea.Field('ensnum',"%f32", "\d*"), ...
                 nmea.Field('pctime', "%f32", "\d*")];
        end
    end
end