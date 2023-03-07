classdef SPDMessage < nmea.Message
% Octans linear speed message
%
%   Defines Octans SPD message with the following fields:
%   x_vel - speed in x direction (m/s)
%   y_vel - speed in y direction (m/s)
%   z_vel - speed in z direction (m/s)
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
        function obj=SPDMessage()
            obj.name = 'SPD';
            obj.msg_id_pattern = obj.name;
            cp = nmea.Field.common_patterns;
            obj.fields = [
                nmea.Field('x_vel', "%f32", cp.float), ...
                nmea.Field('y_vel', "%f32", cp.float),...
                nmea.Field('z_vel', "%f32", cp.float)];
        end
    end
end