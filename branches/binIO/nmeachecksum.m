function isvalid = nmeachecksum(instr)
% isvalid=nmeachecksum(instr)

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

P=inputParser;
P.addRequired('instr',@(x) iscellstr(x) || ischar(x));
P.parse(instr);
infiles=P.Results.instr;

if ischar(infiles)
    instr=cellstr(instr);
end
clear P

%% get position of dollar and asterisk in the string and check if there is a cheksum in the end
tmpdat=textscan([instr{:}],'$ %[^ *] * %2c','Delimiter','');
assert(size(tmpdat{1},1)==size(tmpdat{2},1) && size(tmpdat{2},1)==numel(instr));

%% initialize
checksum = zeros(numel(instr),1,'uint8');
isvalid = false;

%% Calculate checksum
for cntmsg=1:size(tmpdat{1},1)
    NMEA_String_u = uint8(tmpdat{1}{cntmsg});    % convert characters in string to double values
    for count = 1:length(NMEA_String_u)       % checksum calculation ignores $ at start
        checksum(cntmsg) = bitxor(checksum(cntmsg),NMEA_String_u(count));  % checksum calculation
    end
end

%% Check if it matches with checksum in string
% convert checksum to hex value
checksum = dec2hex(checksum,2);

% add leading zero to checksum
isvalid=strcmp(cellstr(checksum),tmpdat{2});
