function [NMEA_GSA,discard]=readGSA(instr)
% readGSA(gsa) interprets a NMEA GSA string
%
%   [GSAstruct]=readGSA(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the gsa data.
%   
%   The output structure will contain the following fields:
%   modesel: Mode selection 0=manual, 1=automatic
%   mode   : Mode 1=No Fix, 2=2D, 3=3D
%   prn    : PRN number of satellites used in solution
%   pdop   : Position dilution of precision
%   hdop   : Horizontal dilution of precision
%   vdop   : Vertical dilution of precision
%
%
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
instr(cellfun(@isempty,instr))=[];
if ischar(instr)
    instr=cellstr(instr);
end

defineNMEA;
%% Scan the data
% format is defined as:
%   $__GSA,a,x,xx,xx,xx,xx,xx,xx,xx,xx,xx,xx,xx,xxx.x,x.x,x.x*hh<CR><LF>
[tmpdat,split]=regexp([instr{:}],patterns.gsa,'names','split');
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
NMEA_GSA.modesel=strcmpi({tmpdat(:).modesel},'A')';
NMEA_GSA.mode=cell2mat(textscan([tmpdat(:).mode],'%u8','delimiter',','));
NMEA_GSA.prn=cell2mat(textscan([tmpdat(:).prn],'%u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8 %u8','delimiter',','));
NMEA_GSA.pdop=cell2mat(textscan([tmpdat(:).pdop],'%f32','delimiter',','));
NMEA_GSA.hdop=cell2mat(textscan([tmpdat(:).hdop],'%f32','delimiter',','));
NMEA_GSA.vdop=cell2mat(textscan([tmpdat(:).vdop],'%f32','delimiter','*'));
