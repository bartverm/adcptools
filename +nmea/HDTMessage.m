classdef HDTMessage < nmea.Message
% HDT True heading NMEA message
%
%   Defines HDT message with the following field:
%   heading - true heading in degrees
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
        function obj=HDTMessage()
            obj.name='HDT';
            obj.msg_id_pattern = obj.name;
            cp = nmea.Field.common_patterns;
            obj.fields = [
                nmea.Field('heading', "%f32", cp.float),...
                nmea.SkipField("T")];
        end
    end
end
