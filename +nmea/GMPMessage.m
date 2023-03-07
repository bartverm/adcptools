classdef GMPMessage < nmea.Message
 % GMP GNSS Map Projection Fix Data NMEA message
 %
 %  GMPMessage defines GMP NMEA message with the following fields:
 %
 %  utc - Nx3 array holding hours, minutes and seconds of time in UTC zone
 %  project_id - either UTM or LOC, the latter being a local projection
 %  project_zone - map projection zone
 %  x - x coordinate
 %  y - y coordinate
 %  mode - nmea.GPSmode array indicating quality of fix. Can have multiple
 %  columns for GPS and GLONASS separately
 %  numsat - Number of satellites in use (00 - 99)
 %  hdop - horizontal dilution of precision using all satellites
 %  alt - antenna altitude with respect to geoid mean sea level (m)
 %  geoid - separation between ellipsoid and geoid (m)
 %  age_dgsp - age of differntial GPS correction (s)
 %  ref_station_id - referece station id
 %
 % see also: nmea, nmea.Message
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
        function obj=GMPMessage()
            obj.name='GMP';
            obj.msg_id_pattern = obj.name;
            cp = nmea.Field.common_patterns;
            obj.fields=[
                nmea.UTCField, ...
                nmea.Field('project_id', "%s", "(?:\w*)?"),...
                nmea.Field('project_zone', "%s", "(?:\w*)?"),...
                nmea.Field('x', "%f64", cp.float),...
                nmea.Field('y', "%f64", cp.float),...
                nmea.Field('mode', "%s", "(?:\w*)?",...
                    @nmea.Field.mode_postproc),...
                nmea.Field('numsat', "%u8", "(?:\d*)?"),...
                nmea.Field('hdop', "%f32", cp.float),...
                nmea.Field('alt', "%f32", cp.float),...
                nmea.Field('geoid', "%f32", cp.float),...
                nmea.Field('age_dgps', "%f32", cp.float),...
                nmea.Field('ref_station_id', "%u16", "(?:\d*)?")];
        end
    end
end