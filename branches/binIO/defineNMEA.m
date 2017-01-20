function [standard, proprietary]=defineNMEA()

% This file defines all regular expressions for nmea strings
% the comf structure contains common patterns encountered in nmea strings
% the patterns structure contains the general nmea string pattern and the
% pattern of several messages
%
% Author: Bart Vermeulen
% Last edit: 02-06-2016

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
comf.float='(?:\-?\d*\.?\d*)'; %optional minus, one or more digits,0 or 1 dot,0 or more digits
comf.time='(?:\d{6}\.?\d*)'; %six digits, 0 or 1 dot, 0 or more digits
comf.csum='[a-fA-F0-9]{2}'; %any of a to f, A to F, 0 to 9, n times
comf.fchr=@(n) ['(?:[a-zA-Z]{',num2str(n),'})']; %any of a to z, A to Z, n times
comf.fdgt=@(n) ['(?:\-?\d{0,',num2str(n),'})']; %optional minus and any of 0-9, n times
comf.status='[AV]'; % A or V
comf.lat='(?:\d{4}\.?\d*)'; %four digits, 0 or 1 dot, 0 or mor digits, N or S
comf.lat_ns='[NS]';
comf.long='(?:\d{5}\.?\d*)'; %five digits, 0 or 1 dot, 0 or mor digits, E or W
comf.long_ew='[EW]';
comf.mag_ew='[EW]';

%% Define standard NMEA patterns
nmea.name='nmea';
nmea.fields(1).name='talker';   nmea.fields(1).pattern=['\$(?<',nmea.fields(1).name,'>',comf.fchr(2),')'];nmea.fields(1).conversion='%s';
nmea.fields(2).name='id';       nmea.fields(2).pattern=['(?<',nmea.fields(2).name,'>',comf.fchr(3),')'];  nmea.fields(2).conversion='%s';
nmea.fields(3).name='content';  nmea.fields(3).pattern=['(?<',nmea.fields(3).name,'>(?:,[^,\*]*)*)\*'];   nmea.fields(3).conversion='%s';
nmea.fields(4).name='checksum'; nmea.fields(4).pattern=['(?<',nmea.fields(4).name,'>',comf.csum,')'];     nmea.fields(4).conversion='%s';
nmea.postaction=[];
standard(1)=nmea;

% global positioning system fix data
gga.name='gga';
gga.fields( 1).name='utc';            gga.fields( 1).pattern=[',(?<',gga.fields( 1).name,'>',comf.time,')?'];         gga.fields( 1).conversion='%2f32 %2f32 %f32'; gga.fields( 1).size=3;    gga.fields( 1).nullval=nan('single');
gga.fields( 2).name='latitude_deg';   gga.fields( 2).pattern=[',(?<',gga.fields( 2).name,'>',comf.lat,')?'];          gga.fields( 2).conversion='%2f64 %f64';       gga.fields( 2).size=2;    gga.fields( 2).nullval=nan('double');
gga.fields( 3).name='lat_ns';         gga.fields( 3).pattern=[',(?<',gga.fields( 3).name,'>',comf.lat_ns,')?'];       gga.fields( 3).conversion='%c';               gga.fields( 3).size=1;    gga.fields( 3).nullval=' ';
gga.fields( 4).name='longitude_deg';  gga.fields( 4).pattern=[',(?<',gga.fields( 4).name,'>',comf.long,')?'];         gga.fields( 4).conversion='%3f64 %f64';       gga.fields( 4).size=2;    gga.fields( 4).nullval=nan('double');
gga.fields( 5).name='lon_ew';         gga.fields( 5).pattern=[',(?<',gga.fields( 5).name,'>',comf.long_ew,')?'];      gga.fields( 5).conversion='%c';               gga.fields( 5).size=1;    gga.fields( 5).nullval=' ';
gga.fields( 6).name='quality';        gga.fields( 6).pattern=[',(?<',gga.fields( 6).name,'>','[0-8]',')?'];           gga.fields( 6).conversion='%u8';              gga.fields( 6).size=1;    gga.fields( 6).nullval=intmax('uint8');
gga.fields( 7).name='num_sat';        gga.fields( 7).pattern=[',(?<',gga.fields( 7).name,'>','(?:[0-1]?[0-9])',')?']; gga.fields( 7).conversion='%u8';              gga.fields( 7).size=1;    gga.fields( 7).nullval=intmax('uint8');
gga.fields( 8).name='hdop';           gga.fields( 8).pattern=[',(?<',gga.fields( 8).name,'>',comf.float,')?'];        gga.fields( 8).conversion='%f32';             gga.fields( 8).size=1;    gga.fields( 8).nullval=nan('single');
gga.fields( 9).name='alt';            gga.fields( 9).pattern=[',(?<',gga.fields( 9).name,'>',comf.float,')?'];        gga.fields( 9).conversion='%f32';             gga.fields( 9).size=1;    gga.fields( 9).nullval=nan('single');
gga.fields(10).name='alt_unit';       gga.fields(10).pattern=[',(?<',gga.fields(10).name,'>',comf.fchr(1),')?'];      gga.fields(10).conversion='%c';               gga.fields(10).size=1;    gga.fields(10).nullval=' ';
gga.fields(11).name='geoid';          gga.fields(11).pattern=[',(?<',gga.fields(11).name,'>',comf.float,')?'];        gga.fields(11).conversion='%f32';             gga.fields(11).size=1;    gga.fields(11).nullval=nan('single');
gga.fields(12).name='geoid_unit';     gga.fields(12).pattern=[',(?<',gga.fields(12).name,'>',comf.fchr(1),')?'];      gga.fields(12).conversion='%c';               gga.fields(12).size=1;    gga.fields(12).nullval=' ';
gga.fields(13).name='age_dgps';       gga.fields(13).pattern=[',(?<',gga.fields(13).name,'>',comf.float,')?'];        gga.fields(13).conversion='%f32';             gga.fields(13).size=1;    gga.fields(13).nullval=nan('single');
gga.fields(14).name='ref_station_id'; gga.fields(14).pattern=[',(?<',gga.fields(14).name,'>',comf.fdgt(4),')?'];      gga.fields(14).conversion='%u16';             gga.fields(14).size=1;    gga.fields(14).nullval=intmax('uint16');
gga.postaction=@post_lat_long;
standard(end+1)=gga;

% depth below transducer
dbt.name='dbt';
dbt.fields(1).name='water_depth_ft'; dbt.fields(1).pattern=[',(?<',dbt.fields(1).name,'>',comf.float  ,')?']; dbt.fields(1).conversion='%f32'; dbt.fields(1).size=1; dbt.fields(1).nullval=nan('single');
dbt.fields(2).name='ft_indicator';   dbt.fields(2).pattern=[',(?<',dbt.fields(2).name,'>',comf.fchr(1),')?']; dbt.fields(2).conversion='%c';   dbt.fields(2).size=1; dbt.fields(2).nullval=(' ');
dbt.fields(3).name='water_depth_m';  dbt.fields(3).pattern=[',(?<',dbt.fields(3).name,'>',comf.float  ,')?']; dbt.fields(3).conversion='%f32'; dbt.fields(3).size=1; dbt.fields(3).nullval=nan('single');
dbt.fields(4).name='m_indicator';    dbt.fields(4).pattern=[',(?<',dbt.fields(4).name,'>',comf.fchr(1),')?']; dbt.fields(4).conversion='%c';   dbt.fields(4).size=1; dbt.fields(4).nullval=(' ');
dbt.fields(5).name='water_depth_f';  dbt.fields(5).pattern=[',(?<',dbt.fields(5).name,'>',comf.float  ,')?']; dbt.fields(5).conversion='%f32'; dbt.fields(5).size=1; dbt.fields(5).nullval=nan('single');
dbt.fields(6).name='f_indicator';    dbt.fields(6).pattern=[',(?<',dbt.fields(6).name,'>',comf.fchr(1),')?']; dbt.fields(6).conversion='%c';   dbt.fields(6).size=1; dbt.fields(6).nullval=(' ');
dbt.postaction=[];
standard(end+1)=dbt;



% dilution of precision (DOP) and active sattelites
gsa.name='gsa';
gsa.fields(1).name='modesel'; gsa.fields(1).pattern=[',(?<',gsa.fields(1).name,'>','[AM]'                       ,')?']; gsa.fields(1).conversion='%c';                                                 gsa.fields(1).size=1;  gsa.fields(1).nullval=nan('single');
gsa.fields(2).name='mode';    gsa.fields(2).pattern=[',(?<',gsa.fields(2).name,'>','[123]'                      ,')?']; gsa.fields(2).conversion='%u8';                                                gsa.fields(2).size=1;  gsa.fields(2).nullval=intmax('uint8');
gsa.fields(3).name='prn';     gsa.fields(3).pattern=[',(?<',gsa.fields(3).name,'>','(?:',comf.fdgt(2),',){12})' ,')?']; gsa.fields(3).conversion='%u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 ';   gsa.fields(3).size=12; gsa.fields(3).nullval=intmax('uint8');
gsa.fields(4).name='pdop';    gsa.fields(4).pattern=[',(?<',gsa.fields(4).name,'>',comf.float                   ,')?']; gsa.fields(4).conversion='%f32';                                               gsa.fields(4).size=1;  gsa.fields(4).nullval=nan('single');
gsa.fields(5).name='hdop';    gsa.fields(5).pattern=[',(?<',gsa.fields(5).name,'>',comf.float                   ,')?']; gsa.fields(5).conversion='%f32';                                               gsa.fields(5).size=1;  gsa.fields(5).nullval=nan('single');
gsa.fields(6).name='vdop';    gsa.fields(6).pattern=[',(?<',gsa.fields(6).name,'>',comf.float                   ,')?']; gsa.fields(6).conversion='%f32';                                               gsa.fields(6).size=1;  gsa.fields(6).nullval=nan('single');
gsa.postaction=[];
standard(end+1)=gsa;

% Track and ground speed
vtg.name='vtg';
vtg.fields(1).name='track_dir_true';        vtg.fields(1).pattern=[',(?<',vtg.fields(1).name,'>',comf.float,')?'  ]; vtg.fields(1).conversion='%f32'; vtg.fields(1).size=1; vtg.fields(1).nullval=nan('single');
vtg.fields(2).name='true_indicator';        vtg.fields(2).pattern=[',(?<',vtg.fields(2).name,'>',comf.fchr(1),')?']; vtg.fields(2).conversion='%c';   vtg.fields(2).size=1; vtg.fields(2).nullval=' ';
vtg.fields(3).name='track_dir_magn';        vtg.fields(3).pattern=[',(?<',vtg.fields(3).name,'>',comf.float,')?'  ]; vtg.fields(3).conversion='%f32'; vtg.fields(3).size=1; vtg.fields(3).nullval=nan('single');
vtg.fields(4).name='magn_indicator';        vtg.fields(4).pattern=[',(?<',vtg.fields(4).name,'>',comf.fchr(1),')?']; vtg.fields(4).conversion='%c';   vtg.fields(4).size=1; vtg.fields(4).nullval=' ';
vtg.fields(5).name='speed_over_ground_kts'; vtg.fields(5).pattern=[',(?<',vtg.fields(5).name,'>',comf.float,')?'  ]; vtg.fields(5).conversion='%f32'; vtg.fields(5).size=1; vtg.fields(5).nullval=nan('single');
vtg.fields(6).name='kts_indicator';         vtg.fields(6).pattern=[',(?<',vtg.fields(6).name,'>',comf.fchr(1),')?']; vtg.fields(6).conversion='%c';   vtg.fields(6).size=1; vtg.fields(6).nullval=' ';
vtg.fields(7).name='speed_over_ground_kmh'; vtg.fields(7).pattern=[',(?<',vtg.fields(7).name,'>',comf.float,')?'  ]; vtg.fields(7).conversion='%f32'; vtg.fields(7).size=1; vtg.fields(7).nullval=nan('single');
vtg.fields(8).name='kmh_indicator';         vtg.fields(8).pattern=[',(?<',vtg.fields(8).name,'>',comf.fchr(1),')?']; vtg.fields(8).conversion='%c';   vtg.fields(8).size=1; vtg.fields(8).nullval=' ';
vtg.fields(9).name='mode_indicator';        vtg.fields(9).pattern=[',(?<',vtg.fields(9).name,'>',comf.fchr(1),')?']; vtg.fields(9).conversion='%c';   vtg.fields(9).size=1; vtg.fields(9).nullval=' ';
vtg.postaction=[];
standard(end+1)=vtg;

% heading
hdt.name='hdt';
hdt.fields(1).name='heading';        hdt.fields(1).pattern=[',(?<',hdt.fields(1).name,'>',comf.float,')?'];   hdt.fields(1).conversion='%f32'; hdt.fields(1).size=1; hdt.fields(1).nullval=nan('single');
hdt.fields(2).name='true_indicator'; hdt.fields(2).pattern=[',(?<',hdt.fields(2).name,'>',comf.fchr(1),')?']; hdt.fields(2).conversion='%c';   hdt.fields(2).size=1; hdt.fields(2).nullval=' ';
hdt.postaction=[];
standard(end+1)=hdt;

% time and date
zda.name='zda';
zda.fields(1).name='utc';  zda.fields(1).pattern=[',(?<',zda.fields(1).name,'>',comf.time,')?'];                                     zda.fields(1).conversion='%2f32 %2f32 %f32';  zda.fields(1).size=3; zda.fields(1).nullval=nan('single');
zda.fields(2).name='date'; zda.fields(2).pattern=[',(?<',zda.fields(2).name,'>',comf.fdgt(2),',',comf.fdgt(2),',',comf.fdgt(4),')']; zda.fields(2).conversion='%2u16,%2u16,%4u16'; zda.fields(2).size=3; zda.fields(2).nullval=intmax('uint16');
zda.fields(3).name='zone'; zda.fields(3).pattern=[',(?<',zda.fields(3).name,'>',comf.fdgt(2),',',comf.fdgt(2),')'];                  zda.fields(3).conversion='%2u8,%2u8';         zda.fields(3).size=2; zda.fields(3).nullval=intmax('uint8');
zda.postaction=[];
standard(end+1)=zda;

% recommended minimum navigation information
rmc.name='rmc';
rmc.fields( 1).name='utc';                     rmc.fields( 1).pattern=[',(?<',rmc.fields( 1).name,'>',comf.time,')?'];            rmc.fields( 1).conversion='%2f32 %2f32 %f32'; rmc.fields( 1).size=3; 	rmc.fields( 1).nullval=nan('single');
rmc.fields( 2).name='status';                  rmc.fields( 2).pattern=[',(?<',rmc.fields( 2).name,'>',comf.status,')?'];          rmc.fields( 2).conversion='%c';               rmc.fields( 2).size=1; 	rmc.fields( 2).nullval=' ';
rmc.fields( 3).name='latitude_deg';            rmc.fields( 3).pattern=[',(?<',rmc.fields( 3).name,'>',comf.lat,')?'];             rmc.fields( 3).conversion='%2f64 %f64';       rmc.fields( 3).size=2;	rmc.fields( 3).nullval=nan('double');
rmc.fields( 4).name='lat_ns';                  rmc.fields( 4).pattern=[',(?<',rmc.fields( 4).name,'>',comf.lat_ns,')?'];          rmc.fields( 4).conversion='%c';               rmc.fields( 4).size=1;	rmc.fields( 4).nullval=' ';
rmc.fields( 5).name='longitude_deg';           rmc.fields( 5).pattern=[',(?<',rmc.fields( 5).name,'>',comf.long,')?'];            rmc.fields( 5).conversion='%2f64 %f64';       rmc.fields( 5).size=2;	rmc.fields( 5).nullval=nan('double');
rmc.fields( 6).name='lon_ew';                  rmc.fields( 6).pattern=[',(?<',rmc.fields( 6).name,'>',comf.long_ew,')?'];         rmc.fields( 6).conversion='%c';               rmc.fields( 6).size=1;	rmc.fields( 6).nullval=' ';
rmc.fields( 7).name='speed_over_ground_kts';   rmc.fields( 7).pattern=[',(?<',rmc.fields( 7).name,'>',comf.float,')?'];           rmc.fields( 7).conversion='%f32';             rmc.fields( 7).size=1; 	rmc.fields( 7).nullval=nan('single');
rmc.fields( 8).name='ground_course';           rmc.fields( 8).pattern=[',(?<',rmc.fields( 8).name,'>',comf.float,')?'];           rmc.fields( 8).conversion='%f32';	            rmc.fields( 8).size=1; 	rmc.fields( 8).nullval=nan('single');
rmc.fields( 9).name='date';                    rmc.fields( 9).pattern=[',(?<',rmc.fields( 9).name,'>',comf.fdgt(6),')?'];         rmc.fields( 9).conversion='%2u8 %2u8 %2u8';   rmc.fields( 9).size=3; 	rmc.fields( 9).nullval=intmax('uint8');
rmc.fields(10).name='magvar';                  rmc.fields(10).pattern=[',(?<',rmc.fields(10).name,'>',comf.float,')?'];           rmc.fields(10).conversion='%f32';	            rmc.fields(10).size=1; 	rmc.fields(10).nullval=nan('single');
rmc.fields(11).name='mag_ew';                  rmc.fields(11).pattern=[',(?<',rmc.fields(11).name,'>',comf.mag_ew,')?'];          rmc.fields(11).conversion='%c';               rmc.fields(11).size=1;	rmc.fields(11).nullval=' ';
rmc.fields(12).name='mode';                    rmc.fields(12).pattern=[',?(?<',rmc.fields(12).name,'>[ADEMSN]?)?'];               rmc.fields(12).conversion='%c';               rmc.fields(12).size=1;	rmc.fields(12).nullval=' '; % in some description it seems the status is not always output, so made optional
rmc.postaction=[];
standard(end+1)=rmc;



% % geographic position
% gll.name='gll';
% gll.pattern=['\$\w{2}GLL,',...
%               '(?<lat>',comf.lat,',)',...
%               '(?<long>',comf.long,',)'...
%               '(?<utc>',comf.time,',)',...
%               '(?<status>',comf.status,')[\*,]',...
%               '(?<mode>[ADEMSN]?)?\*?',...
%               comf.fhex(2)];
% patterns(end+1)=gll;




% % depth below surface
% dbs.name='dbs';
% dbs.pattern=['\$\w{2}DBS,',...
%               '(?<depthf>',comf.float,',)',...
%               'f,',...
%               '(?<depthM>',comf.float,',)',...
%               'M,',...
%               '(?<depthF>',comf.float,',)',...
%               'F\*',...
%               comf.fhex(2)];
% patterns(end+1)=dbs;




%% Define proprietary NMEA strings

% RDI ensemble number and pc time
rdens.name='rdens';
rdens.head='\$RDENS';
rdens.fields(1).name='ensnum'; rdens.fields(1).pattern=[',(?<',rdens.fields(1).name,'>\d*)?']; rdens.fields(1).conversion='%u32'; rdens.fields(1).size=1; rdens.fields(1).nullval=intmax('uint32');
rdens.fields(2).name='pctime'; rdens.fields(2).pattern=[',(?<',rdens.fields(2).name,'>\d*)?']; rdens.fields(2).conversion='%u32'; rdens.fields(2).size=1; rdens.fields(2).nullval=intmax('uint32');
rdens.tail=',PC';
rdens.tail_is_checksum=false;
rdens.postaction=[];
proprietary(1)=rdens;

% % heading, pitch and roll
% hpr.name='hpr';
% hpr.pattern=['\$PSAT,HPR,'...
%               '(?<utc>',comf.time,',)',...
%               '(?<heading>',comf.float,',)',...
%               '(?<pitch>',comf.float,',)',...
%               '(?<roll>',comf.float,',)',...
%               '\w\*',...
%               comf.fhex(2)];
% patterns(end+1)=hpr;

% % sattelite fault detection
% gbs.name='gbs';
% gbs.pattern=['\$PSAT,GBS,',...
%               '(?<utc>',comf.time,',)',...
%               '(?<errlat>',comf.float,',)',...
%               '(?<errlong>',comf.float,',)',...
%               '(?<erralt>',comf.float,',)',...
%               '(?<failsat>',comf.fdgt(2),',)',...
%               '(?<phprfail>',comf.float,',)',...
%               '(?<rngbias>',comf.float,',)',...
%               '(?<sdrngbias>',comf.float,',)',...
%               '(?<flag>[012]?\*)',...
%               comf.fhex(2)];
% patterns(end+1)=gbs;
end



function out=post_lat_long(in)
       out=in;
       out.latitude=out.latitude_deg(1,:)+out.latitude_deg(2,:)/60;
       out.latitude(out.lat_ns=='S')=-out.latitude(out.lat_ns=='S');
       out.longitude=out.longitude_deg(1,:)+out.longitude_deg(2,:)/60;
       out.longitude(out.lon_ew=='W')=-out.longitude(out.lon_ew=='W');       
end



