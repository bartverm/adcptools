function [NMEA_ZDA,discard]=readZDA(instr)
% readZDA(zda) interprets a NMEA ZDA string
%
%   [ZDAstruct]=readZDA(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the zda data.
%   
%   The output structure will contain the following fields:
%   UTCtime : Hour minutes and seconds of the UTC time
%   date    : Date in year month day
%   zone    : Timezone expressed in timeshift in hours and minutes
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


% tmpdat=textscan([instr{:}],'$ %*2s ZDA %2f32 %2f32 %f32 %u16 %u16 %u16 %d8 %d8 %*2s',nlines,'Delimiter',',*');
defineNMEA;
[tmpdat,split]=regexp([instr{:}],patterns.zda,'names','split');
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
NMEA_GGA.UTCtime=cell2mat(textscan([tmpdat(:).utc],'%2f32 %2f32 %f32','delimiter',','));
NMEA_ZDA.date=fliplr(cell2mat(textscan([tmpdat(:).date],'%u16 %u16 %u16','delimiter',',')));
NMEA_ZDA.zone=cell2mat(textscan([tmpdat(:).zone],'%d8 %d8','delimiter',',*'));
