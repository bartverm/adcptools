classdef TROMessage < nmea.Message
% Octans tilt rotation message
%
%   Defines TRO Octans message with the following fields:
%   pitch - in degrees, negative for bow down
%   roll - in degrees negative for port down
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
        function obj=TROMessage()
            obj.name = 'TRO';
            obj.msg_id_pattern = obj.name;
            cp = nmea.Field.common_patterns;
            obj.fields = [
                nmea.Field('pitch',["%f32", "%s"],...
                    "(?:" + cp.float + ",\w)?" ,...
                    @(x) nmea.TROMessage.tilt_postproc(x, 'P')),...
                nmea.Field('roll',["%f32", "%s"],...
                    "(?:" + cp.float + ",\w)?" ,...
                    @(x) nmea.TROMessage.tilt_postproc(x, 'B'))];
        end
    end
    methods (Static, Access=protected)
        function out=tilt_postproc(in,neg_char)
            out=in{1};
            fneg=strcmp(in{2},neg_char);
            out(fneg)=-out(fneg);
        end
    end
end