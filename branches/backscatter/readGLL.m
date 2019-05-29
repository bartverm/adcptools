function [NMEA_GLL,discard]=readGLL(instr)
% readGLL(gll) interprets a NMEA GLL string
%
%   [GLLstruct]=readGLL(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the gll data.
%   
%   The output structure will contain the following fields:
%   UTCtime : Hour minutes and seconds of the UTC time
%   lat     : Latitude in decimal degrees. Negative for the
%             southern hemisphere
%   long    : Longitude in decimal degrees. Negative for eastern longitude
%   status  : true=valid, false=invalid
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
[tmpdat,split]=regexp([instr{:}],patterns.gll,'names','split');
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

%Initialize variable
NMEA_GLL.UTCtime=cell2mat(textscan([tmpdat(:).utc],'%2f32 %2f32 %f32','delimiter',','));
tmplat=textscan([tmpdat(:).lat],'%2f64 %f64 %s','delimiter',',');
NMEA_GLL.lat=tmplat{1}+tmplat{2}/60;
isneg=strcmp(tmplat{3},'S');
NMEA_GLL.lat(isneg)=-NMEA_GLL.lat(isneg);
clear tmplat
tmplong=textscan([tmpdat(:).long],'%3f64 %f64 %s','delimiter',',');
NMEA_GLL.long=tmplong{1}+tmplong{2}/60;
isneg=strcmp(tmplong{3},'W');
NMEA_GLL.long(isneg)=-NMEA_GLL.long(isneg);
clear tmplong isneg
NMEA_GLL.status=strcmp({tmpdat(:).status},'A')';
if isfield(tmpdat,'mode')
    NMEA_GLL.mode=zeros(size(tmpdat),'uint8')';
    tmpmod={tmpdat(:).mode};
    NMEA_GLL.mode(strcmp(tmpmod,'A'))=1;
    NMEA_GLL.mode(strcmp(tmpmod,'D'))=2;
    NMEA_GLL.mode(strcmp(tmpmod,'E'))=3;
    NMEA_GLL.mode(strcmp(tmpmod,'M'))=4;
    NMEA_GLL.mode(strcmp(tmpmod,'S'))=5;
end
