classdef VTGMessage < nmea.Message
% VTG Course over ground and ground speed data
%
%   Defines VTG message with the following fields:
%   track_dir_true - Course over ground in degrees True
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
        function obj=VTGMessage()
            obj.name='VTG';
            obj.msg_id_pattern = obj.name;
            cp = nmea.Field.common_patterns;
            obj.fields=[
                nmea.Field('track_dir_true', "%f32", cp.float),...
                nmea.SkipField("T?"),...
                nmea.Field('track_dir_magn', "%f32", cp.float),...
                nmea.SkipField("M?"),...
                nmea.Field('speed_over_ground_kts', "%f32", cp.float),...
                nmea.SkipField("N?"),...
                nmea.Field('speed_over_ground_kmh', "%f32", cp.float),...
                nmea.SkipField("K?"),...
                nmea.Field('mode_indicator', "%s", "(?:\w*)?",...
                    @nmea.Field.mode_postproc)];
        end
    end
end
