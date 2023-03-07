classdef GPSMode < uint8
% Specifies the quality or mode of acquired GPS data
%
%   GPSMode members:
%   Unknown
%   NoFix 
%   Autonomous
%   Differential
%   Precise
%   RealTimeKinematic
%   FloatRTK
%   Estimated
%   Manual
%   Simulated
%
% see also: nmea, nmea.Message
%
%   Copyright 2020, Bart Vermeulen

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


    enumeration
        Unknown (255)
        NoFix (78) %N
        Autonomous (65) %A
        Differential (68) %D
        Precise (80) %P
        RealTimeKinematic (82) %R
        FloatRTK (70) %F
        Estimated (69) %E
        Manual (77) %M
        Simulated (83) %S
    end
end