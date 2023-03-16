classdef LINMessage < nmea.Message
% Octans linear motion data
%
%   Octans LIN message with the following fields:
%   surge_x - surge motion in m
%   sway_y - sway motion in m
%   heave_z - heave motion in m
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
        function obj=LINMessage()
            obj.name='LIN';
            obj.msg_id_pattern = obj.name;
            cp = nmea.Field.common_patterns;
            obj.fields=[nmea.Field('surge_x',"%f32", cp.float), ...
                 nmea.Field('sway_y', "%f32", cp.float),...
                 nmea.Field('heave_z', "%f32", cp.float)];
        end
    end
end