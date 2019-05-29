function [NMEA_GBS,discard]=readGBS(instr)
% readGBS(gbs) interprets the csi propietary NMEA message GBS
%
%   [GBSstruct]=readGBS(instr) Reads a NMEA sentence in the character array
%   or cell of strings and returns a structure with the gbs data.
%   
%   The output structure will contain the following fields:
%   UTCtime : Hour minutes and seconds of the UTC time
%   errlat  : Estimated error in latitude
%   errlong : Estimated error in longitude
%   erralt  : Estimated error in altitude
%   failID  : ID of satallite that most likely failed
%   pHPRfail: Probability of HPR failure
%   rbias   : Estimate of range bias on most likely failed satellite
%   biasstd : Standard deviation of bias estimate
%   flag    : Good (0) / Warning (1) / Bad (2) Flag
%
%   Warning: This script was not tested yet...
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

[tmpdat,split]=regexp([instr{:}],patterns.gbs,'names','split');
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

% tmpdat=textscan([instr{:}],'$PSAT GBS %2f32 %2f32 %f32 %f32 %f32 %f32 %u16 %f32 %f32 %f32 %u8 %*2s',nlines,'Delimiter',',*');



%Initialize variable
NMEA_GBS.UTCtime=cell2mat(textscan([tmpdat(:).utc],'%2f32 %2f32 %f32','delimiter',','));
NMEA_GBS.errlat=cell2mat(textscan([tmpdat(:).errlat],'%f32','delimiter',','));
NMEA_GBS.errlong=cell2mat(textscan([tmpdat(:).errlong],'%f32','delimiter',','));
NMEA_GBS.erralt=cell2mat(textscan([tmpdat(:).erralt],'%f32','delimiter',','));
NMEA_GBS.failID=cell2mat(textscan([tmpdat(:).failsat],'%u16','delimiter',','));
NMEA_GBS.pHPRfail=cell2mat(textscan([tmpdat(:).phprfail],'%f32','delimiter',','));
NMEA_GBS.rbias=cell2mat(textscan([tmpdat(:).rngbias],'%f32','delimiter',','));
NMEA_GBS.biasstd=cell2mat(textscan([tmpdat(:).sdrngbias],'%f32','delimiter',','));
NMEA_GBS.flag=cell2mat(textscan([tmpdat(:).flag],'%u8','delimiter','*'));
