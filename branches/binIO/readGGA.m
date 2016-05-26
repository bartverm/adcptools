function [NMEA_GGA,discard]=readGGA(instr)
% readGGA(gga) interprets a NMEA GGA string
%
%   [GGAstruct]=readGGA(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the gga data.
%   
%   The output structure will contain the following fields:
%   UTCtime : Hour minutes and seconds of the UTC time
%   lat     : Latitude in decimal degrees. Negative for the
%             southern hemisphere
%   long    : Longitude in decimal degrees. Negative for western longitude
%   qualind : A quality indicator defined as follows:
%                 0 = fix not available or invalid
%                 1 = GPS fix
%                 2 = Differential GPS fix
%                 3 = GPS PPS Mode, fix valid
%                 4 = Real Time Kinematic. System used in RTK mode with 
%                     fixed integers
%                 5 = Float RTK. Satellite system used in RTK mode, 
%                     floating integers
%                 6 = Estimated (dead reckoning) mode
%                 7 = Manual Input Mode
%                 8 = Simulator mode
%    numsat : Number of satellites in use
%    hdop   : Horizontal Dilution of Precisition
%    antalt : Antenna altitude in meters
%
%    Author: Bart Vermeulen
%    Date last edit: 21-12-2009

%    Copyright 2009 Bart Vermeulen
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

% Create input parser
P=inputParser;
P.addRequired('instr',@(x) (iscellstr(x) | ischar(x)));
P.parse(instr);

% Convert input to cell of strings if it is a character array
if ischar(instr)
    instr=cellstr(instr);
end

defineNMEA

% tmpdat=textscan([instr{:}],'$ %*2s GGA %2f32 %2f32 %f32 %2f64 %f64 %s %3f64 %f64 %s %u8 %u8 %f32 %f32 %*1s %f32 %*1s %f32 %u16 %*2s',nlines,'Delimiter',',*');
[tmpdat,split]=regexp([instr{:}],patterns.gga,'names','split');
clear instr
discard=find(~strcmp(split,''));
fdiscard=[];
if any(discard)
    for cdis=1:length(discard)
        nbad=numel(regexp(split{discard(cdis)},patterns.nmea));
        fdiscard=[fdiscard,repmat(discard(cdis),1,nbad)]; %#ok<AGROW>
    end
    discard=fdiscard+(1:length(fdiscard))-1;     
end

%read the data
NMEA_GGA.UTCtime=cell2mat(textscan([tmpdat(:).utc],'%2f32 %2f32 %f32','delimiter',','));
tmplat=textscan([tmpdat(:).lat],'%2f64 %f64 %s','delimiter',',');
NMEA_GGA.lat=tmplat{1}+tmplat{2}/60;
isneg=strcmp(tmplat{3},'S');
NMEA_GGA.lat(isneg)=-NMEA_GGA.lat(isneg);
clear tmplat
tmplong=textscan([tmpdat(:).long],'%3f64 %f64 %s','delimiter',',');
NMEA_GGA.long=tmplong{1}+tmplong{2}/60;
isneg=strcmp(tmplong{3},'W');
NMEA_GGA.long(isneg)=-NMEA_GGA.long(isneg);
clear tmplong isneg
NMEA_GGA.qualind=cell2mat(textscan([tmpdat(:).qualind],'%u8','delimiter',','));
NMEA_GGA.numsat=cell2mat(textscan([tmpdat(:).numsat],'%u8','delimiter',','));
NMEA_GGA.hdop=cell2mat(textscan([tmpdat(:).hdop],'%f32','delimiter',','));
NMEA_GGA.antalt=cell2mat(textscan([tmpdat(:).antalt],'%f32','delimiter',','));
NMEA_GGA.geosep=cell2mat(textscan([tmpdat(:).geosep],'%f32','delimiter',','));
NMEA_GGA.agediff=cell2mat(textscan([tmpdat(:).agediff],'%f32','delimiter',','));
NMEA_GGA.diffid=cell2mat(textscan([tmpdat(:).diffid],'%u16','delimiter','*'));


