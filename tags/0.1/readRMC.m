function [NMEA_RMC,discard]=readRMC(instr)
% readRMC(rmc) interprets a NMEA RMC string
%
%   [RMCstruct]=readRMC(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the rmc data.
%   
%   The output structure will contain the following fields:
%   UTCtime : Hour minutes and seconds of the UTC time
%   valid   : 0 invalid, 1 valid
%   lat     : Latitude in decimal degrees. Negative for the
%             southern hemisphere
%   long    : Longitude in decimal degrees. Negative for western longitude
%   grspeed : Ground speed in knots
%   trmgood : Track made good referenced to the north
%   magvar  : Magnetic variation in degrees. Negative for western variation
%   date    : Year month and day
%   mode    : Indicates the mode: 1 autonomous, 2 differential, 3
%             estimated, 4 manual, 5 simulator, 0 not valid
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

defineNMEA;

[tmpdat,split]=regexp([instr{:}],patterns.rmc,'names','split');
% clear instr
discard=find(~strcmp(split,''));
fdiscard=[];
if any(discard)
    for cdis=1:length(discard)
        nbad=numel(regexp(split{discard(cdis)},patterns.nmea));
        fdiscard=[fdiscard,repmat(discard(cdis),1,nbad)]; %#ok<AGROW>
    end
    discard=fdiscard+(1:length(fdiscard))-1;     
end



% tmpdat=textscan([instr{:}],'$ %*2s RMC %2f32 %2f32 %f32 %1s %2f64 %f64 %s %3f64 %f64 %s %f32 %f32 %2u8 %2u8 %2u8 %f32 %1s %1s %*2s',nlines,'Delimiter',',*');
%Initialize variable
NMEA_RMC.UTCtime=cell2mat(textscan([tmpdat(:).utc],'%2f32 %2f32 %f32','delimiter',','));
NMEA_RMC.status=strcmp({tmpdat(:).status},'A')';
tmplat=textscan([tmpdat(:).lat],'%2f64 %f64 %s','delimiter',',');
NMEA_RMC.lat=tmplat{1}+tmplat{2}/60;
isneg=strcmp(tmplat{3},'S');
NMEA_RMC.lat(isneg)=-NMEA_RMC.lat(isneg);
clear tmplat
tmplong=textscan([tmpdat(:).long],'%3f64 %f64 %s','delimiter',',');
NMEA_RMC.long=tmplong{1}+tmplong{2}/60;
isneg=strcmp(tmplong{3},'W');
NMEA_RMC.long(isneg)=-NMEA_RMC.long(isneg);
clear tmplong isneg
NMEA_RMC.grspeed=cell2mat(textscan([tmpdat(:).groundSpKn],'%f32','delimiter',','));
NMEA_RMC.trmgood=cell2mat(textscan([tmpdat(:).groundCourse],'%f32','delimiter',','));
NMEA_RMC.date=fliplr(cell2mat(textscan([tmpdat(:).date],'%2u8 %2u8 %2u8','delimiter',',')));
tmpmv=textscan([tmpdat(:).magvar],'%f64 %s','delimiter',',*');
NMEA_RMC.magvar=tmpmv{1};
isneg=strcmp(tmpmv{2},'W');
NMEA_RMC.magvar(isneg)=-NMEA_RMC.magvar(isneg);
clear tmpmv isneg
if isfield(tmpdat,'mode')
    NMEA_RMC.mode=zeros(size(tmpdat),'uint8')';
    tmpmod={tmpdat(:).mode};
    NMEA_RMC.mode(strcmp(tmpmod,'A'))=1;
    NMEA_RMC.mode(strcmp(tmpmod,'D'))=2;
    NMEA_RMC.mode(strcmp(tmpmod,'E'))=3;
    NMEA_RMC.mode(strcmp(tmpmod,'M'))=4;
    NMEA_RMC.mode(strcmp(tmpmod,'S'))=5;
end

