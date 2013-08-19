function [NMEA_VTG,discard]=readVTG(instr)
% readVTG(vtg) interprets a NMEA VTG string
%
%   [VTGstruct]=readVTG(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the vtg data.
%   
%   The output structure will contain the following fields:
%   TrackDegTrue : Track, degrees true
%   TrackDegMagn : Track, degrees magentic
%   SpeedKnots   : Speed in knots
%   SpeedKmH     : Speed in Km/h
%
%    Author: Bart Vermeulen
%    Date last edit: 17-12-2009

%    Copyright 2009,2010 Bart Vermeulen, Frans Buschman
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

%    21-09-2010: patterns.hdt changed to patterns.vtg (Frans)

% Create input parser
P=inputParser;
P.addRequired('instr',@(x) (iscellstr(x) | ischar(x)));
P.parse(instr);

% Convert input to cell of strings if it is a character array
if ischar(instr)
    instr=cellstr(instr);
end

% tmpdat=textscan([instr{:}],'$ %*2s VTG %f32 T %f32 M %f32 N %f32 K %*2s','Delimiter',nlines,',*');
defineNMEA;
[tmpdat,split]=regexp([instr{:}],patterns.vtg,'names','split'); 
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
NMEA_VTG.TrackDegTrue=cell2mat(textscan([tmpdat(:).TrueCourse],'%f32','delimiter',','));
NMEA_VTG.TrackDegMagn=cell2mat(textscan([tmpdat(:).MagCourse],'%f32','delimiter',','));
NMEA_VTG.SpeedKnots=cell2mat(textscan([tmpdat(:).speedKn],'%f32','delimiter',','));
NMEA_VTG.SpeedKmH=cell2mat(textscan([tmpdat(:).speedKm],'%f32','delimiter',','));
if isfield(tmpdat,'mode')
    NMEA_VTG.mode=zeros(size(tmpdat),'uint8')';
    tmpmod={tmpdat(:).mode};
    NMEA_VTG.mode(strcmp(tmpmod,'A'))=1;
    NMEA_VTG.mode(strcmp(tmpmod,'D'))=2;
    NMEA_VTG.mode(strcmp(tmpmod,'E'))=3;
    NMEA_VTG.mode(strcmp(tmpmod,'M'))=4;
    NMEA_VTG.mode(strcmp(tmpmod,'S'))=5;
end
