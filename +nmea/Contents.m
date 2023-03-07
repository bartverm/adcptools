% %%%%%%%%%%%%
% NMEA toolbox
% %%%%%%%%%%%%
%
% Copyright 2020 Bart Vermeulen
%
% This package contains classes to process NMEA strings. To get started,
% this is how it works:
%
% 1. Load the available nmea messages:
%       msgs=nmea.Message.get_all_messages();
% 
% 2. Parse string or character array data:
%       dat=msgs.parse(str);
%
% You can find more detailed explanation below:
%
% NMEA classes:
%   nmea.Message - Defines messages and parses NMEA data
%   nmea.Field - Defines data fields in NMEA messages
%   nmea.GPSMode - Enumeration specifying the kind of GPS fix
%
% If you want to add support for new messages, you have to subclass the
% nmea.Message class. The easiest way to start is using one of the existing
% NMEA messages (e.g. nmea.GGAMessage) copy it and modify it for the NMEA
% message you want to read. Also, once the class works, add it to the
% Message.get_all_messages function



%     This file is part of the NMEA toolbox.
% 
%     NMEA toolbox is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     NMEA toolbox is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with NMEA toolbox.  If not, see <https://www.gnu.org/licenses/>.
