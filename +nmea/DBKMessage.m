classdef DBKMessage < nmea.DBTMessage
% depth below keel message
    methods
        function obj = DBKMessage()
            obj = obj@nmea.DBTMessage;
            obj.name = 'DBK';
            obj.msg_id_pattern = obj.name;
        end
    end
end