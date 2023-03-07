classdef ZDAMessage < nmea.Message
% ZDA Time and date NMEA message
%
%   Defines ZDA message with the following field:
%   datetime - datetime object array holding time and timezone
%
%   see also: nmea, nmea.Message
%
%   Copyright 2020, Bart Vermeulen

%     This file is part of the NMEA toolbox.
% 
%     NMEA toolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     NMEA toolbox is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with NMEA toolbox.  If not, see <https://www.gnu.org/licenses/>.

    methods
        function obj = ZDAMessage()
            obj.name='ZDA';
            obj.msg_id_pattern = obj.name;
            cp = nmea.Field.common_patterns;
            obj.fields = nmea.Field('datetime',...
                ["%2f32" "%2f32" "%f32" "%2u16" "%2u16",...
                "%4u16" "%2d8" "%2u8"],...
                cp.time + "," + ...
                cp.ndgt(2) + "," + cp.ndgt(2) + "," + cp.ndgt(4) + "," + ...
                cp.ndgt(2) + "," + cp.ndgt(2),...
                @nmea.ZDAMessage.datetime_postproc);
        end
    end

    methods (Static, Access=protected)
        function out=datetime_postproc(in)
            tzone=unique(reshape(sprintf('%+03d%02u',...
                in{7},in{8}),5,[])','rows');
            if size(tzone,1)>1
                warning(['Timezone changing in ZDA data, using:',...
                    tzone(1,:)])
            end
            out=datetime(in{6}, in{5}, in{4}, in{1}, in{2}, in{3},...
                'Timezone',tzone(1,:));
        end
    end
end

