classdef GGAMessage < nmea.Message
% GGA Global Positioning System Fix Data
%
% Defines nmea GGA message with the following fields:
%   utc - Nx3 array with hours, minutes and seconds UTC time
%   latitude - in decimal degrees, negative for south
%   longitude - in decimal degrees, negative for west
%   quality - nmea.GPSMode array indicating quality of GPS data
%   numsat - number of satellites in use (0 - 12)
%   hdop - horizontal dilution of precision
%   alt - altitude wrt geoid msl (m)
%   geoid - separation between WGS84 ellipsoid and geoid msl (m)
%   age_dgps - Time in seconds since last DGPS update
%   ref_station_id - Differential reference station ID (0000-1023)
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
        function obj = GGAMessage()
            obj.msg_id_pattern = 'GGA';
            obj.name = 'GGA';
            cp = nmea.Field.common_patterns;
            obj.fields=[...
                nmea.UTCField, ...
                nmea.LatitudeField,...
                nmea.LongitudeField,...
                nmea.Field('quality',"%u8","[0-8]",...
                    @nmea.GGAMessage.quality_postproc),...
                nmea.Field('numsat', "%u8", "(?:[0-9]?[0-9]?)?"),...
                nmea.Field('hdop', "%f32", cp.float),...
                nmea.Field('alt', "%f32", cp.float),...
                nmea.SkipField("M?"),...
                nmea.Field('geoid', "%f32", cp.float),...
                nmea.SkipField("M?"),...
                nmea.Field('age_dgps', "%f32", cp.float),...
                nmea.Field('ref_station_id',"%u16", cp.ndgt(4))];
        end
    end
    methods (Static, Access=protected)
        function out=quality_postproc(in)
            out(numel(in{1}),1)=nmea.GPSMode.Unknown;
            out(in{1}==0)=nmea.GPSMode.NoFix;
            out(in{1}==1)=nmea.GPSMode.Autonomous;
            out(in{1}==2)=nmea.GPSMode.Differential;
            out(in{1}==3)=nmea.GPSMode.Precise;
            out(in{1}==4)=nmea.GPSMode.RealTimeKinematic;
            out(in{1}==5)=nmea.GPSMode.FloatRTK;
            out(in{1}==6)=nmea.GPSMode.Estimated;
            out(in{1}==7)=nmea.GPSMode.Manual;
            out(in{1}==8)=nmea.GPSMode.Simulated;
        end
    end
end