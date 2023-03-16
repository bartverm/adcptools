classdef INFMessage < nmea.Message
% Octans status message
%
%   INF message holding the field:
%   octans_status - holds the octans status in hexadecimal byte values
%
%   see also: nmea, nmea.Message
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
        function obj=INFMessage()
            obj.name='INF';
            obj.msg_id_pattern = obj.name;
            obj.fields=nmea.Field('octans_status', "%s", "(?:\w*)?");
        end
    end
end