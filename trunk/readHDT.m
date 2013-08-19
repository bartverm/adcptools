function [NMEA_HDT,discard]=readHDT(instr)
% readHDT(hdt) interprets a NMEA HDT string
%
%   [HDTstruct]=readHDT(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the hdt data.
%   
%   The output structure will contain the following fields:
%   heading: heading of the gps compass
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
[tmpdat,split]=regexp([instr{:}],patterns.hdt,'names','split');
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
%tmpdat=textscan([instr{:}],'$ %*2s HDT %f32 T %*2s',nlines,'Delimiter',',*');

%Initialize variable
NMEA_HDT.heading=cell2mat(textscan([tmpdat(:).heading],'%f32','delimiter',','));

