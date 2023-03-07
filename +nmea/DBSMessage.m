classdef DBSMessage < nmea.DBTMessage
% depth below sea surface message
    methods
        function obj = DBSMessage()
            obj = obj@nmea.DBTMessage;
            obj.name = 'DBS';
            obj.msg_id_pattern = obj.name;
        end
    end
end