% This file defines all regular expressions for nmea strings
% the comf structure contains common patterns encountered in nmea strings
% the patterns structure contains the general nmea string pattern and the
% pattern of several messages
%
% Author: Bart Vermeulen
% Last edit: 13-01-2010

%    Copyright 2010 Bart Vermeulen
%
%    This file is part of ADCPTools.
%
%    ADCPTools is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    ADCPTools is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with ADCPTools.  If not, see <http://www.gnu.org/licenses/>.


%% Define common fields (they are all optional, so empty field is allowed. Except the checksum) 
comf.float='(?:\-?\d*\.?\d*)?'; %optional minus, one or more digits,0 or 1 dot,0 or more digits 
comf.time='(?:\d{6}\.?\d*)?'; %six digits, 0 or 1 dot, 0 or more digits
comf.fhex=@(n) ['[a-fA-F0-9]{',num2str(n),'}']; %any of a to f, A to F, 0 to 9, n times
comf.fchr=@(n) ['(?:[a-zA-Z]{',num2str(n),'})?']; %any of a to z, A to Z, n times
comf.fdgt=@(n) ['(?:\-?\d{0,',num2str(n),'})?']; %optional minus and any of 0-9, n times
comf.status='[AV]?'; % A or V
comf.lat='(?:\d{4}\.?\d*)?,[NS]?'; %four digits, 0 or 1 dot, 0 or mor digits, N or S
comf.long='(?:\d{5}\.?\d*)?,[EW]?'; %four digits, 0 or 1 dot, 0 or mor digits, E or W

%% Define patterns
patterns.nmea='\$(\w{2})(\w*),(\w+)[^\*\f\n\r\t\v]+(?:\*\w{2}|PC)';
% global positioning system fix data
patterns.gga=['\$\w{2}GGA,'...
              '(?<utc>',comf.time,',)',...
              '(?<lat>',comf.lat,',)',...
              '(?<long>',comf.long,',)'...
              '(?<qualind>[0-8],)'...
              '(?<numsat>(?:[0-9]?[0-9])?,)',...
              '(?<hdop>',comf.float,',)',...
              '(?<antalt>',comf.float,',)',...
              'M?,',...
              '(?<geosep>',comf.float,',)',...
              'M?,',...
              '(?<agediff>',comf.float,',)',...
              '(?<diffid>',comf.fdgt(4),'\*)',...
              comf.fhex(2)];
% dilution of precision (DOP) and active sattelites
patterns.gsa=['\$\w{2}GSA,',...
              '(?<modesel>[AM]?),',...
              '(?<mode>[123]?,)',...
              '(?<prn>(?:',comf.fdgt(2),',){12})',...
              '(?<pdop>',comf.float,',)',...
              '(?<hdop>',comf.float,',)',...
              '(?<vdop>',comf.float,'\*)',...
              comf.fhex(2)];
% RDI ensemble number and pc time
patterns.rdens=['\$RDENS,',...
                '(?<ensnum>\d*,)',...
                '(?<pctime>\d*,)',...
                'PC'];
% sattelite fault detection
patterns.gbs=['\$PSAT,GBS,',...
              '(?<utc>',comf.time,',)',...
              '(?<errlat>',comf.float,',)',...
              '(?<errlong>',comf.float,',)',...
              '(?<erralt>',comf.float,',)',...
              '(?<failsat>',comf.fdgt(2),',)',...
              '(?<phprfail>',comf.float,',)',...
              '(?<rngbias>',comf.float,',)',...
              '(?<sdrngbias>',comf.float,',)',...
              '(?<flag>[012]?\*)',...
              comf.fhex(2)];
% geographic position
patterns.gll=['\$\w{2}GLL,',...
              '(?<lat>',comf.lat,',)',...
              '(?<long>',comf.long,',)'...
              '(?<utc>',comf.time,',)',...
              '(?<status>',comf.status,')[\*,]',...
              '(?<mode>[ADEMSN]?)?\*?',...
              comf.fhex(2)];
% depth below transducer
patterns.dbt=['\$\w{2}DBT,',...
              '(?<depthf>',comf.float,',)',...
              'f,',...
              '(?<depthM>',comf.float,',)',...
              'M,',...
              '(?<depthF>',comf.float,',)',...
              'F\*',...
              comf.fhex(2)];    
% heading, pitch and roll
patterns.hpr=['\$PSAT,HPR,'...
              '(?<utc>',comf.time,',)',...
              '(?<heading>',comf.float,',)',...
              '(?<pitch>',comf.float,',)',...
              '(?<roll>',comf.float,',)',...
              '\w\*',...
              comf.fhex(2)];
% depth below surface
patterns.dbs=['\$\w{2}DBS,',...
              '(?<depthf>',comf.float,',)',...
              'f,',...
              '(?<depthM>',comf.float,',)',...
              'M,',...
              '(?<depthF>',comf.float,',)',...
              'F\*',...
              comf.fhex(2)]; 
% time and date
patterns.zda=['\$\w{2}ZDA,',...
              '(?<utc>',comf.time,',)',...
              '(?<date>',comf.fdgt(2),',',comf.fdgt(2),',',comf.fdgt(4),',)',...
              '(?<zone>',comf.fdgt(2),',',comf.fdgt(2),'\*)',...
              comf.fhex(2)];    
% recommended minimum navigation information
patterns.rmc=['\$\w{2}RMC,',...
              '(?<utc>',comf.time,',)',...
              '(?<status>',comf.status,'),',...
              '(?<lat>',comf.lat,',)',...
              '(?<long>',comf.long,',)'...
              '(?<groundSpKn>',comf.float,',)',...
              '(?<groundCourse>',comf.float,',)',...
              '(?<date>',comf.fdgt(6),',)',...
              '(?<magvar>',comf.float,',[EW]?[,\*])',...
              '(?<mode>[ADEMSN]?)?\*?',...
              comf.fhex(2)];
% heading
patterns.hdt=['\$\w{2}HDT,',...
              '(?<heading>',comf.float,',)',...
              'T\*',...
              comf.fhex(2)];
% Track and ground speed
patterns.vtg=['\$\w{2}VTG,'...
              '(?<TrueCourse>',comf.float,',)',...
              'T,',...
              '(?<MagCourse>',comf.float,',)',...
              'M?,',...
              '(?<speedKn>',comf.float,',)',...
              'N?,',...
              '(?<speedKm>',comf.float,',)',...
              'K[,\*]',...
              '(?<mode>[ADEMSN]?)?\*?',...
              comf.fhex(2)];

